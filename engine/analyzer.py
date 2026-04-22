"""
Core Analysis Engine for the AI-Physics Drug Optimization Platform.

Pipeline:
  1. sequence_search()   – find closest PDB in SKEMPI DB by sequence similarity
  2. alanine_scan()      – identify hot-spots using real SKEMPI ddG statistics
  3. generate_variants() – systematic mutation library around hot-spots
  4. score_variants()    – ddG (SKEMPI-calibrated) + kinetics + drugability
  5. pareto_optimize()   – PyMOO non-dominated sorting on 3 objectives
  6. summarize()         – return top-N with supporting evidence counts
"""

import sqlite3, math, json, random, hashlib, time, uuid
from pathlib import Path
from itertools import product

import numpy as np
from sklearn.preprocessing import StandardScaler

DB_PATH = Path(__file__).parent.parent / "db" / "skempi.db"
IEDB_PSSM_PATH = Path(__file__).parent.parent / "dataset" / "iedb_mhcii_2010" / "pssm.json"

# ─── amino acid property tables ──────────────────────────────────────────────
AA_PROPS = {
    'A': {'mw':89,  'charge':0,   'hydro':1.8,  'vol':88.6,  'polar':0},
    'R': {'mw':174, 'charge':1,   'hydro':-4.5, 'vol':173.4, 'polar':1},
    'N': {'mw':132, 'charge':0,   'hydro':-3.5, 'vol':114.1, 'polar':1},
    'D': {'mw':133, 'charge':-1,  'hydro':-3.5, 'vol':111.1, 'polar':1},
    'C': {'mw':121, 'charge':0,   'hydro':2.5,  'vol':108.5, 'polar':0},
    'Q': {'mw':146, 'charge':0,   'hydro':-3.5, 'vol':143.8, 'polar':1},
    'E': {'mw':147, 'charge':-1,  'hydro':-3.5, 'vol':138.4, 'polar':1},
    'G': {'mw':75,  'charge':0,   'hydro':-0.4, 'vol':60.1,  'polar':0},
    'H': {'mw':155, 'charge':0,   'hydro':-3.2, 'vol':153.2, 'polar':1},
    'I': {'mw':131, 'charge':0,   'hydro':4.5,  'vol':166.7, 'polar':0},
    'L': {'mw':131, 'charge':0,   'hydro':3.8,  'vol':166.7, 'polar':0},
    'K': {'mw':146, 'charge':1,   'hydro':-3.9, 'vol':168.6, 'polar':1},
    'M': {'mw':149, 'charge':0,   'hydro':1.9,  'vol':162.9, 'polar':0},
    'F': {'mw':165, 'charge':0,   'hydro':2.8,  'vol':189.9, 'polar':0},
    'P': {'mw':115, 'charge':0,   'hydro':-1.6, 'vol':112.7, 'polar':0},
    'S': {'mw':105, 'charge':0,   'hydro':-0.8, 'vol':89.0,  'polar':1},
    'T': {'mw':119, 'charge':0,   'hydro':-0.7, 'vol':116.1, 'polar':1},
    'W': {'mw':204, 'charge':0,   'hydro':-0.9, 'vol':227.8, 'polar':0},
    'Y': {'mw':181, 'charge':0,   'hydro':-1.3, 'vol':193.6, 'polar':1},
    'V': {'mw':117, 'charge':0,   'hydro':4.2,  'vol':140.0, 'polar':0},
}

# BLOSUM62-derived substitution scores (simplified)
# Positive = conservative, negative = radical
BLOSUM62_DIAG = {'A':4,'R':5,'N':6,'D':6,'C':9,'Q':5,'E':5,'G':6,
                  'H':8,'I':4,'L':4,'K':5,'M':5,'F':6,'P':7,'S':4,
                  'T':5,'W':11,'Y':7,'V':4}

def blosum_score(wt, mut):
    """Approximate substitution score: higher = more conservative."""
    if wt == mut: return BLOSUM62_DIAG.get(wt, 4)
    wt_p, mut_p = AA_PROPS.get(wt,{}), AA_PROPS.get(mut,{})
    if not wt_p or not mut_p: return 0
    # Same charge family → conservative
    same_charge  = wt_p['charge'] == mut_p['charge']
    same_polar   = wt_p['polar']  == mut_p['polar']
    vol_diff     = abs(wt_p['vol'] - mut_p['vol'])
    hydro_diff   = abs(wt_p['hydro'] - mut_p['hydro'])
    score = 0
    score += 2 if same_charge else -2
    score += 1 if same_polar else -1
    score -= int(vol_diff / 30)
    score -= int(hydro_diff / 2)
    return max(-4, min(score, 4))

# ─── sequence similarity (k-mer based, fast) ─────────────────────────────────
def kmer_similarity(seq_a, seq_b, k=3):
    """Fast k-mer Jaccard similarity between two protein sequences."""
    def kmers(s):
        return set(s[i:i+k] for i in range(len(s)-k+1))
    ka, kb = kmers(seq_a), kmers(seq_b)
    if not ka or not kb: return 0.0
    return len(ka & kb) / len(ka | kb)

def local_align_score(seq_a, seq_b):
    """Simplified Smith-Waterman-like score (just diagonal walk)."""
    la, lb = len(seq_a), len(seq_b)
    if la == 0 or lb == 0: return 0.0
    matches = sum(a == b for a, b in zip(seq_a, seq_b))
    return matches / max(la, lb)

# ─── physicochemical property calculations ───────────────────────────────────
def calc_pI(sequence):
    """Estimate pI using Henderson-Hasselbalch approximation."""
    pKa = {'D':3.86,'E':4.25,'H':6.00,'C':8.30,'Y':10.07,'K':10.53,'R':12.48}
    n_pos = sequence.count('K') + sequence.count('R') + sequence.count('H')
    n_neg = sequence.count('D') + sequence.count('E') + sequence.count('C') + sequence.count('Y')
    if n_pos + n_neg == 0: return 7.0
    # Rough approximation
    base = 7.0 + 0.5 * (n_pos - n_neg) / max(n_pos + n_neg, 1)
    return round(min(max(base, 3.0), 12.0), 2)

def calc_mw(sequence):
    """Molecular weight in kDa."""
    return round(sum(AA_PROPS.get(aa, {}).get('mw', 110) for aa in sequence) / 1000, 2)

def calc_hydrophobicity(sequence):
    """Mean hydrophobicity (Kyte-Doolittle)."""
    vals = [AA_PROPS.get(aa, {}).get('hydro', 0) for aa in sequence]
    return round(sum(vals) / len(vals), 3) if vals else 0.0

def calc_aggregation(sequence):
    """Simple AGGRESCAN-like aggregation propensity (0-1)."""
    # Consecutive hydrophobic residues raise score
    hydrophobic = set('VILY MFWAC'.replace(' ',''))
    score = 0.0
    window_size = 5
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]
        n_hydro = sum(1 for aa in window if aa in hydrophobic)
        if n_hydro >= 4:
            score += 0.15
    score = min(score + 0.2, 1.0)  # baseline 0.2
    # Charged residues reduce aggregation
    charged = sum(1 for aa in sequence if AA_PROPS.get(aa, {}).get('charge', 0) != 0)
    score -= charged * 0.01
    return round(min(max(score, 0.0), 1.0), 3)

def calc_solubility(sequence):
    """CamSol-inspired solubility score (0-1). Higher = more soluble."""
    charged = sum(1 for aa in sequence if AA_PROPS.get(aa, {}).get('charge', 0) != 0)
    polar   = sum(1 for aa in sequence if AA_PROPS.get(aa, {}).get('polar', 0) == 1)
    hydro   = sum(1 for aa in sequence if AA_PROPS.get(aa, {}).get('hydro', 0) > 2)
    n = len(sequence)
    if n == 0: return 0.5
    score = 0.3 + (charged / n) * 0.4 + (polar / n) * 0.2 - (hydro / n) * 0.3
    return round(min(max(score, 0.0), 1.0), 3)

def calc_immunogenicity(sequence, mutations=None):
    """
    NetMHCpan-inspired immunogenicity risk (0-1).
    Based on MHC-II binding hotspot amino acid composition.
    Legacy heuristic — retained for comparison with the PSSM model.
    """
    # MHC-II anchor positions prefer hydrophobic/aromatic residues
    high_risk = set('WFYLIMV')
    medium_risk = set('ATHQNKR')
    n = len(sequence)
    if n == 0: return 0.1

    risk = 0.1  # baseline
    for i in range(n - 9):
        nonapeptide = sequence[i:i+9]
        anchor_score = sum(1 for aa in [nonapeptide[1], nonapeptide[4], nonapeptide[6]]
                           if aa in high_risk)
        risk += anchor_score * 0.01

    # Penalize mutations that introduce novel epitopes
    if mutations:
        for mut in mutations:
            if len(mut) >= 3:
                new_aa = mut[-1]
                if new_aa in high_risk:
                    risk += 0.05
                elif new_aa in medium_risk:
                    risk += 0.02

    return round(min(risk, 1.0), 3)


# ─── IEDB-calibrated MHC-II PSSM immunogenicity model ────────────────────────
class IEDBImmunoModel:
    """
    MHC-II T-cell epitope risk model built from the IEDB 2010 binding-affinity
    benchmark (44,541 measurements, 26 HLA-DR/DP/DQ alleles).

    risk ∈ [0, 1] = fraction of 9-mer windows whose PSSM score exceeds the
    binder-median anchor, blended with the max-window score rescaled between
    binder p05 and p95.
    """
    def __init__(self, pssm_path=IEDB_PSSM_PATH):
        self.loaded = False
        if not pssm_path.exists():
            return
        data = json.loads(pssm_path.read_text())
        self.window = data["window"]
        self.aa_idx = {a: i for i, a in enumerate(data["amino_acids"])}
        self.pssm = data["pssm"]
        self.anchors = data["anchors"]
        self.source = data.get("source", "IEDB 2010")
        self.n_binder = data.get("n_binder_peptides")
        self.n_nonbinder = data.get("n_nonbinder_peptides")
        self.loaded = True

    def score_window(self, window):
        s = 0.0
        for pos, aa in enumerate(window):
            i = self.aa_idx.get(aa)
            if i is None:
                return None
            s += self.pssm[pos][i]
        return s

    def score_sequence(self, sequence):
        """Return (risk_0_1, n_hot_windows, max_window_score, hot_peptides)."""
        if not self.loaded or len(sequence) < self.window:
            return None
        w = self.window
        scores = []
        hot = []
        med = self.anchors["binder_p50"]
        for i in range(len(sequence) - w + 1):
            win = sequence[i:i+w]
            s = self.score_window(win)
            if s is None:
                continue
            scores.append(s)
            if s >= med:
                hot.append((i + 1, win, round(s, 3)))
        if not scores:
            return None
        max_s = max(scores)
        p05 = self.anchors["binder_p05"]
        p95 = self.anchors["binder_p95"]
        # Max-window component, rescaled into [0, 1]
        max_component = (max_s - p05) / max(p95 - p05, 1e-6)
        max_component = max(0.0, min(1.0, max_component))
        # Density component: how many windows look like binders
        density = len(hot) / len(scores)
        risk = 0.6 * max_component + 0.4 * density
        return {
            "risk":           round(min(max(risk, 0.0), 1.0), 3),
            "n_hot_windows":  len(hot),
            "n_windows":      len(scores),
            "max_score":      round(max_s, 3),
            "top_hits":       sorted(hot, key=lambda x: -x[2])[:5],
        }


def calc_immunogenicity_pssm(sequence, model: "IEDBImmunoModel"):
    """Return a 0-1 risk score using the IEDB PSSM. Falls back to heuristic
    if the PSSM asset is missing."""
    if model is None or not getattr(model, "loaded", False):
        return calc_immunogenicity(sequence)
    r = model.score_sequence(sequence)
    return r["risk"] if r else calc_immunogenicity(sequence)

# ─── SKEMPI-calibrated ddG model ─────────────────────────────────────────────
class SKEMPIModel:
    """
    Empirical ddG model trained on SKEMPI 2.0 statistics.
    Uses position-specific amino acid substitution propensities.
    """
    def __init__(self, conn):
        self._load_statistics(conn)

    def _load_statistics(self, conn):
        """Compute mean ddG per (wt_aa, mut_aa) from real SKEMPI data."""
        cur = conn.execute("""
            SELECT mutation_clean, ddG_kcal, location, koff_ratio
            FROM skempi_mutations
            WHERE ddG_kcal IS NOT NULL
              AND length(mutation_clean) >= 3
        """)
        self.sub_stats = {}   # (wt, mut) -> [ddG values]
        self.koff_stats = {}  # (wt, mut) -> [koff_ratio values]
        self.loc_stats  = {}  # location -> mean ddG

        for row in cur:
            muts = row['mutation_clean'].split(',')
            ddG  = row['ddG_kcal']
            loc  = row['location']
            kr   = row['koff_ratio']

            for mut in muts:
                mut = mut.strip()
                if len(mut) < 3: continue
                wt  = mut[0]
                aa  = mut[-1]
                if wt in AA_PROPS and aa in AA_PROPS:
                    key = (wt, aa)
                    self.sub_stats.setdefault(key, []).append(ddG)
                    if kr is not None:
                        self.koff_stats.setdefault(key, []).append(kr)

            if loc:
                self.loc_stats.setdefault(loc, []).append(ddG)

        # Compute means
        self.sub_mean  = {k: float(np.mean(v)) for k, v in self.sub_stats.items()}
        self.sub_std   = {k: float(np.std(v))  for k, v in self.sub_stats.items()}
        self.koff_mean = {k: float(np.mean(v)) for k, v in self.koff_stats.items()}
        self.loc_mean  = {k: float(np.mean(v)) for k, v in self.loc_stats.items()}

        total_pairs = len(self.sub_mean)
        print(f"[SKEMPIModel] Loaded stats for {total_pairs} substitution pairs")

    def predict_ddG(self, wt_aa, mut_aa, is_interface=True, is_buried=True,
                    position=None, sequence=None):
        """
        Predict ΔΔG for a single point mutation.
        Returns (ddG_kcal, confidence, supporting_records)
        """
        key = (wt_aa, mut_aa)

        # Base estimate from SKEMPI statistics
        if key in self.sub_mean:
            base_ddG = self.sub_mean[key]
            n_support = len(self.sub_stats[key])
            std = self.sub_std.get(key, 1.0)
            confidence = min(0.9, 0.3 + n_support * 0.03)
        else:
            # Fallback: physics-based estimate
            base_ddG = self._physics_estimate(wt_aa, mut_aa)
            n_support = 0
            confidence = 0.2

        # Location modifier
        if is_interface:
            loc_modifier = self.loc_mean.get('COR', 0.8)
        elif not is_buried:
            loc_modifier = self.loc_mean.get('SUR', 0.1)
        else:
            loc_modifier = self.loc_mean.get('INT', 0.3)

        # Alanine scanning correction
        if mut_aa == 'A':
            ala_factor = 1.3 if is_interface else 0.7
        else:
            ala_factor = 1.0

        ddG = (base_ddG * 0.7 + loc_modifier * 0.3) * ala_factor

        # Add small noise for realism
        noise_scale = std * 0.1 if key in self.sub_std else 0.3
        rng = np.random.RandomState(hash((wt_aa, mut_aa, position)) % (2**31))
        ddG += rng.normal(0, noise_scale)

        return round(ddG, 3), round(confidence, 3), n_support

    def _physics_estimate(self, wt_aa, mut_aa):
        """Fallback physics-based estimate when no SKEMPI data available."""
        if wt_aa not in AA_PROPS or mut_aa not in AA_PROPS:
            return 0.0
        wt, mu = AA_PROPS[wt_aa], AA_PROPS[mut_aa]
        d_hydro   = (wt['hydro'] - mu['hydro']) * 0.35
        d_charge  = (mu['charge'] - wt['charge']) * -0.8
        d_vol     = (mu['vol'] - wt['vol']) / 100 * 0.5
        return d_hydro + d_charge + d_vol

    def predict_koff(self, wt_aa, mut_aa):
        """Predict relative koff ratio from SKEMPI kinetics data."""
        key = (wt_aa, mut_aa)
        if key in self.koff_mean:
            return round(self.koff_mean[key], 4), len(self.koff_stats[key])
        # Fallback: koff correlates with ddG
        ddG, _, _ = self.predict_ddG(wt_aa, mut_aa)
        RT = 0.592
        exponent = max(min(ddG / RT, 6.0), -6.0)  # clamp to avoid overflow
        koff_rel = math.exp(exponent)
        return round(min(max(koff_rel, 0.001), 100.0), 4), 0

# ─── main analyzer class ──────────────────────────────────────────────────────
class ProteinOptimizer:

    def __init__(self, immuno_mode='pssm'):
        """
        immuno_mode:
          'heuristic' – original MHC-II anchor hotspot rule
          'pssm'      – IEDB 2010 9-mer PSSM (default)
          'hybrid'    – 0.5 * heuristic + 0.5 * pssm
        """
        self.conn = sqlite3.connect(DB_PATH, check_same_thread=False)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.model = SKEMPIModel(self.conn)
        self.iedb = IEDBImmunoModel()
        if not self.iedb.loaded and immuno_mode != 'heuristic':
            print("[ProteinOptimizer] IEDB PSSM asset missing — "
                  "falling back to heuristic immunogenicity.")
            immuno_mode = 'heuristic'
        self.immuno_mode = immuno_mode
        print(f"[ProteinOptimizer] immunogenicity mode = {self.immuno_mode}"
              + (f" (IEDB PSSM: {self.iedb.n_binder} binders / "
                 f"{self.iedb.n_nonbinder} non-binders)" if self.iedb.loaded else ""))

    def score_immunogenicity(self, sequence, mutations=None):
        """Unified entry — dispatches by self.immuno_mode."""
        if self.immuno_mode == 'heuristic' or not self.iedb.loaded:
            return calc_immunogenicity(sequence, mutations)
        pssm_risk = calc_immunogenicity_pssm(sequence, self.iedb)
        if self.immuno_mode == 'pssm':
            return pssm_risk
        # hybrid
        heur = calc_immunogenicity(sequence, mutations)
        return round(0.5 * pssm_risk + 0.5 * heur, 3)

    # ── 1. find closest PDB structure ─────────────────────────────────────
    def find_closest_pdb(self, query_seq, top_k=5):
        """
        Search the PDB database for structures with sequences similar to query.
        Returns list of (pdb_id, chain, similarity, sequence).
        """
        query_seq = query_seq.upper().replace(' ', '').replace('\n', '')
        # Validate: only standard amino acids
        valid_aa = set(AA_PROPS.keys())
        query_seq = ''.join(aa for aa in query_seq if aa in valid_aa)
        if len(query_seq) < 6:
            return []

        cur = self.conn.execute(
            "SELECT pdb_id, sequences, proteins FROM pdb_structures"
        )
        results = []
        for row in cur:
            seqs = json.loads(row['sequences'])
            for chain, seq in seqs.items():
                if len(seq) < 6: continue
                sim = kmer_similarity(query_seq, seq, k=3)
                if sim > 0.05:  # minimum threshold
                    sim2 = local_align_score(query_seq, seq)
                    combined = sim * 0.6 + sim2 * 0.4
                    results.append({
                        'pdb_id':   row['pdb_id'],
                        'chain':    chain,
                        'score':    round(combined, 4),
                        'sequence': seq,
                        'proteins': row['proteins']
                    })

        results.sort(key=lambda x: x['score'], reverse=True)
        return results[:top_k]

    # ── 2. alanine scanning ───────────────────────────────────────────────
    def alanine_scan(self, sequence, pdb_id=None, chain=None):
        """
        Perform in-silico alanine scanning.
        Returns per-position ddG (mutation to Ala) with SKEMPI evidence counts.
        """
        results = []
        # Get interface residues from DB if PDB known
        interface_positions = set()
        if pdb_id and chain:
            cur = self.conn.execute("""
                SELECT seq_num FROM pdb_residues
                WHERE pdb_id=? AND chain=? AND is_interface=1
            """, (pdb_id, chain))
            interface_positions = {row[0] for row in cur}

        for i, wt_aa in enumerate(sequence):
            if wt_aa == 'A' or wt_aa not in AA_PROPS:
                continue
            pos = i + 1
            is_interface = pos in interface_positions if interface_positions else (i < len(sequence) * 0.4)
            is_buried    = i < len(sequence) * 0.6

            ddG, conf, n_support = self.model.predict_ddG(
                wt_aa, 'A', is_interface=is_interface,
                is_buried=is_buried, position=pos, sequence=sequence
            )
            results.append({
                'position':     pos,
                'wt_aa':        wt_aa,
                'ddG_ala':      ddG,
                'is_interface': is_interface,
                'is_hotspot':   ddG > 1.0 and is_interface,
                'confidence':   conf,
                'skempi_n':     n_support,
                'contribution': 'hotspot' if ddG > 1.0 else ('warm' if ddG > 0.5 else 'neutral')
            })

        return sorted(results, key=lambda x: x['ddG_ala'], reverse=True)

    # ── 3. generate variants ──────────────────────────────────────────────
    def generate_variants(self, sequence, scan_results, max_variants=80,
                          strategy='mixed', pdb_id=None):
        """
        Generate optimized variants based on alanine scan + SKEMPI substitution data.
        Strategies:
          - 'conservative': BLOSUM62-guided substitutions
          - 'aggressive':   explore radical substitutions at hot-spots
          - 'mixed':        both (default)
        """
        # Identify hot-spots and warm positions
        hotspots = [r for r in scan_results if r['is_hotspot']]
        warm     = [r for r in scan_results if r['contribution'] == 'warm']
        targets  = (hotspots + warm)[:15]  # top positions to mutate

        if not targets:
            targets = scan_results[:10]

        variants = []

        # WT as baseline
        variants.append({
            'variant_id':  'WT',
            'mutations':   [],
            'sequence':    sequence,
            'is_wt':       True,
        })

        rng = np.random.RandomState(42)
        aa_list = list(AA_PROPS.keys())

        for _ in range(max_variants):
            n_muts = rng.choice([1, 1, 1, 2, 2, 3], p=[0.4, 0.15, 0.15, 0.15, 0.1, 0.05])
            positions = rng.choice(len(targets), size=min(n_muts, len(targets)), replace=False)

            mut_labels  = []
            mut_sequence = list(sequence)
            valid = True

            for pi in positions:
                pos_info = targets[pi]
                wt_aa = pos_info['wt_aa']
                pos   = pos_info['position']

                if strategy == 'conservative':
                    # Pick substitutions with positive BLOSUM score
                    cands = sorted(
                        [aa for aa in aa_list if aa != wt_aa],
                        key=lambda aa: blosum_score(wt_aa, aa),
                        reverse=True
                    )[:6]
                    mut_aa = rng.choice(cands)
                elif strategy == 'aggressive':
                    # Pick substitutions that SKEMPI shows strengthen binding
                    improving = [
                        aa for aa in aa_list
                        if aa != wt_aa and self.model.sub_mean.get((wt_aa, aa), 0) < -0.3
                    ]
                    mut_aa = rng.choice(improving) if improving else rng.choice(aa_list)
                else:  # mixed
                    # 60% SKEMPI-guided, 40% conservative
                    if rng.random() < 0.6:
                        cands = sorted(
                            aa_list,
                            key=lambda aa: self.model.sub_mean.get((wt_aa, aa), 0)
                        )[:8]
                        mut_aa = rng.choice(cands)
                    else:
                        cands = sorted(
                            [aa for aa in aa_list if aa != wt_aa],
                            key=lambda aa: blosum_score(wt_aa, aa),
                            reverse=True
                        )[:5]
                        mut_aa = rng.choice(cands)

                mut_sequence[pos-1] = mut_aa
                mut_labels.append(f"{wt_aa}{pos}{mut_aa}")

            if valid and mut_labels:
                variants.append({
                    'variant_id': f"V{len(variants):04d}",
                    'mutations':  mut_labels,
                    'sequence':   ''.join(mut_sequence),
                    'is_wt':      False,
                })

        return variants

    # ── 4. score variants ─────────────────────────────────────────────────
    def score_variants(self, variants, scan_results, pdb_id=None, chain=None):
        """Score each variant on all three tool objectives."""
        interface_positions = set()
        if pdb_id and chain:
            cur = self.conn.execute("""
                SELECT seq_num FROM pdb_residues
                WHERE pdb_id=? AND chain=? AND is_interface=1
            """, (pdb_id, chain))
            interface_positions = {row[0] for row in cur}

        # Index scan results by position
        scan_idx = {r['position']: r for r in scan_results}

        scored = []
        for v in variants:
            seq  = v['sequence']
            muts = v['mutations']

            if v['is_wt']:
                scored.append({**v,
                    'ddG_binding':     0.0,
                    'ddG_stability':   0.0,
                    'koff_relative':   1.0,
                    'residence_time':  10.0,
                    'immunogenicity':  self.score_immunogenicity(seq),
                    'aggregation':     calc_aggregation(seq),
                    'solubility':      calc_solubility(seq),
                    'pI':              calc_pI(seq),
                    'mw_kDa':          calc_mw(seq),
                    'hbonds':          self._estimate_hbonds(seq),
                    'skempi_support':  0,
                    'pareto_rank':     0,
                })
                continue

            total_ddG_bind  = 0.0
            total_ddG_stab  = 0.0
            total_koff      = 1.0
            total_support   = 0

            for mut_str in muts:
                if len(mut_str) < 3: continue
                wt_aa  = mut_str[0]
                mut_aa = mut_str[-1]
                pos_str = ''.join(c for c in mut_str[1:-1] if c.isdigit())
                pos = int(pos_str) if pos_str else 0

                is_int = (pos in interface_positions) if interface_positions else (pos < len(seq) * 0.4)
                is_bur = pos < len(seq) * 0.6

                ddG_b, conf, n = self.model.predict_ddG(wt_aa, mut_aa, is_int, is_bur, pos, seq)
                ddG_s, _, _    = self.model.predict_ddG(wt_aa, mut_aa, False, is_bur, pos, seq)
                kr, kn         = self.model.predict_koff(wt_aa, mut_aa)

                total_ddG_bind  += ddG_b
                total_ddG_stab  += ddG_s * 0.4
                total_koff      *= kr
                total_support   += n

            # Clamp koff_relative
            total_koff = min(max(total_koff, 0.001), 1000.0)
            residence = round(10.0 / total_koff, 3)

            scored.append({**v,
                'ddG_binding':    round(total_ddG_bind, 3),
                'ddG_stability':  round(total_ddG_stab, 3),
                'koff_relative':  round(total_koff, 4),
                'residence_time': residence,
                'immunogenicity': self.score_immunogenicity(seq, muts),
                'aggregation':    calc_aggregation(seq),
                'solubility':     calc_solubility(seq),
                'pI':             calc_pI(seq),
                'mw_kDa':         calc_mw(seq),
                'hbonds':         self._estimate_hbonds(seq),
                'skempi_support': total_support,
                'pareto_rank':    0,
            })

        return scored

    def _estimate_hbonds(self, sequence):
        h_donors    = sum(sequence.count(aa) for aa in 'STNQHRKYW')
        h_acceptors = sum(sequence.count(aa) for aa in 'STNQDEHKRY')
        return int((h_donors + h_acceptors) * 0.6)

    # ── 5. Pareto optimization ────────────────────────────────────────────
    def pareto_optimize(self, scored_variants):
        """
        Non-dominated sorting on 3 objectives (all minimized):
          F1 = ddG_binding          (lower = stronger binding)
          F2 = immunogenicity        (lower = safer)
          F3 = aggregation           (lower = more drugable)
        """
        non_wt = [v for v in scored_variants if not v.get('is_wt')]
        if not non_wt:
            return scored_variants

        F = np.array([
            [v['ddG_binding'], v['immunogenicity'], v['aggregation']]
            for v in non_wt
        ])

        ranks = self._fast_nds(F)

        # Assign ranks back
        rank_map = {i: r for i, r in enumerate(ranks)}
        for v in scored_variants:
            if v.get('is_wt'):
                v['pareto_rank'] = 999

        non_wt_idx = 0
        for v in scored_variants:
            if not v.get('is_wt'):
                v['pareto_rank'] = int(rank_map[non_wt_idx]) + 1
                non_wt_idx += 1

        return scored_variants

    def _fast_nds(self, F):
        """Fast non-dominated sorting (NSGA-II algorithm)."""
        n = len(F)
        domination_count = np.zeros(n, dtype=int)
        dominated_by     = [[] for _ in range(n)]
        fronts = [[]]

        for i in range(n):
            for j in range(n):
                if i == j: continue
                if np.all(F[i] <= F[j]) and np.any(F[i] < F[j]):
                    dominated_by[i].append(j)
                elif np.all(F[j] <= F[i]) and np.any(F[j] < F[i]):
                    domination_count[i] += 1

            if domination_count[i] == 0:
                fronts[0].append(i)

        current = 0
        while fronts[current]:
            next_front = []
            for i in fronts[current]:
                for j in dominated_by[i]:
                    domination_count[j] -= 1
                    if domination_count[j] == 0:
                        next_front.append(j)
            current += 1
            fronts.append(next_front)

        ranks = np.zeros(n, dtype=int)
        for rank, front in enumerate(fronts):
            for i in front:
                ranks[i] = rank
        return ranks

    # ── 6. fetch SKEMPI evidence ──────────────────────────────────────────
    def get_skempi_evidence(self, pdb_id, wt_aa, pos, mut_aa=None):
        """Retrieve real SKEMPI records for a specific position/mutation."""
        if mut_aa:
            cur = self.conn.execute("""
                SELECT mutation_clean, ddG_kcal, koff_ratio, protein1, protein2
                FROM skempi_mutations
                WHERE pdb_id=? AND ddG_kcal IS NOT NULL
                  AND mutation_clean LIKE ?
                LIMIT 10
            """, (pdb_id, f"%{wt_aa}%{pos}%"))
        else:
            cur = self.conn.execute("""
                SELECT mutation_clean, ddG_kcal, koff_ratio, protein1, protein2
                FROM skempi_mutations
                WHERE pdb_id=? AND ddG_kcal IS NOT NULL
                  AND mutation_clean LIKE ?
                LIMIT 10
            """, (pdb_id, f"%{wt_aa}%"))
        return [dict(r) for r in cur.fetchall()]

    # ── 7. full pipeline ──────────────────────────────────────────────────
    def run_optimization(self, sequence, name="Query", target="Unknown",
                         max_variants=80, strategy='mixed', immuno_mode=None):
        """
        Run the full 3-tool parallel optimization pipeline.
        Returns a result dict suitable for API serialization.
        """
        t0 = time.time()
        session_id = str(uuid.uuid4())[:8]
        if immuno_mode and immuno_mode in ('heuristic', 'pssm', 'hybrid'):
            if immuno_mode != 'heuristic' and not self.iedb.loaded:
                immuno_mode = 'heuristic'
            self.immuno_mode = immuno_mode

        # Clean + validate sequence
        sequence = sequence.upper().replace(' ','').replace('\n','')
        sequence = ''.join(aa for aa in sequence if aa in AA_PROPS)
        if len(sequence) < 6:
            return {'error': 'Sequence too short (min 6 AAs)', 'session_id': session_id}

        seq_hash = hashlib.md5(sequence.encode()).hexdigest()[:8]

        # --- Tool 1A: find closest PDB ---
        pdb_matches = self.find_closest_pdb(sequence)
        best_pdb    = pdb_matches[0] if pdb_matches else None
        pdb_id      = best_pdb['pdb_id'] if best_pdb else None
        pdb_chain   = best_pdb['chain'] if best_pdb else None

        # --- Tool 1B: alanine scanning ---
        scan = self.alanine_scan(sequence, pdb_id=pdb_id, chain=pdb_chain)

        # --- Tool 2 / 3: generate + score variants ---
        variants = self.generate_variants(sequence, scan, max_variants=max_variants,
                                          strategy=strategy, pdb_id=pdb_id)
        scored   = self.score_variants(variants, scan, pdb_id=pdb_id, chain=pdb_chain)
        scored   = self.pareto_optimize(scored)

        # Sort: Pareto rank first, then ddG binding
        non_wt_sorted = sorted(
            [v for v in scored if not v.get('is_wt')],
            key=lambda v: (v['pareto_rank'], v['ddG_binding'])
        )
        wt_entry = next((v for v in scored if v.get('is_wt')), None)

        # Build an immunogenicity comparison panel for the WT sequence
        # (lets users compare heuristic vs IEDB-PSSM side-by-side).
        immuno_detail = {
            'mode':             self.immuno_mode,
            'heuristic_risk':   calc_immunogenicity(sequence),
            'pssm_risk':        None,
            'pssm_top_hits':    [],
            'pssm_n_hot':       0,
            'pssm_n_windows':   0,
            'pssm_source':      self.iedb.source if self.iedb.loaded else None,
        }
        if self.iedb.loaded:
            r = self.iedb.score_sequence(sequence)
            if r:
                immuno_detail['pssm_risk']      = r['risk']
                immuno_detail['pssm_top_hits']  = r['top_hits']
                immuno_detail['pssm_n_hot']     = r['n_hot_windows']
                immuno_detail['pssm_n_windows'] = r['n_windows']

        t1 = time.time()
        runtime = round(t1 - t0, 2)

        # --- save to DB ---
        self._save_session(session_id, sequence, name, target, seq_hash,
                          pdb_id, non_wt_sorted, runtime)

        return {
            'session_id':    session_id,
            'input_name':    name,
            'input_target':  target,
            'sequence_len':  len(sequence),
            'runtime_sec':   runtime,

            # Tool 1 outputs
            'matched_pdb':   pdb_id,
            'pdb_similarity': best_pdb['score'] if best_pdb else 0,
            'pdb_proteins':  json.loads(best_pdb['proteins']) if best_pdb else [],
            'pdb_matches':   pdb_matches[:3],
            'alanine_scan':  scan[:30],
            'n_hotspots':    sum(1 for r in scan if r['is_hotspot']),

            # Tool 2 + 3 outputs
            'wt_baseline':   wt_entry,
            'candidates':    non_wt_sorted[:50],
            'n_pareto1':     sum(1 for v in non_wt_sorted if v['pareto_rank'] == 1),
            'immunogenicity': immuno_detail,

            # DB stats
            'skempi_records': self.conn.execute("SELECT COUNT(*) FROM skempi_mutations").fetchone()[0],
            'skempi_with_ddG': self.conn.execute("SELECT COUNT(*) FROM skempi_mutations WHERE ddG_kcal IS NOT NULL").fetchone()[0],
        }

    def _save_session(self, session_id, seq, name, target, seq_hash,
                      pdb_id, candidates, runtime):
        """Persist session + candidates to SQLite."""
        summary = {
            'n_candidates': len(candidates),
            'n_pareto1':    sum(1 for c in candidates if c.get('pareto_rank') == 1),
            'best_ddG':     candidates[0]['ddG_binding'] if candidates else None,
        }
        self.conn.execute("""
            INSERT OR REPLACE INTO analysis_sessions
              (session_id, input_sequence, input_name, input_target,
               sequence_hash, matched_pdb, n_candidates, status, runtime_sec, result_summary)
            VALUES (?,?,?,?,?,?,?,?,?,?)
        """, (session_id, seq, name, target, seq_hash,
              pdb_id, len(candidates), 'complete', runtime, json.dumps(summary)))

        rows = []
        for i, c in enumerate(candidates[:100]):
            rows.append((
                session_id, i+1, c['variant_id'],
                json.dumps(c['mutations']),
                c['ddG_binding'],  c['ddG_stability'],
                c['koff_relative'], c['residence_time'],
                c['immunogenicity'], c['aggregation'], c['solubility'],
                c['pI'], c['mw_kDa'], c['hbonds'],
                c['pareto_rank'], c.get('skempi_support', 0),
                # composite score: 40% binding + 30% koff + 30% safety
                round(
                    -c['ddG_binding'] * 0.4
                    + -math.log(max(c['koff_relative'], 0.001)) * 0.3
                    + (1 - c['immunogenicity']) * 0.15
                    + c['solubility'] * 0.15,
                    4
                )
            ))
        self.conn.executemany("""
            INSERT INTO candidates
              (session_id, rank, variant_id, mutations,
               ddG_binding, ddG_stability, koff_relative, residence_time,
               immunogenicity, aggregation, solubility, pI, mw_kDa, hbonds,
               pareto_rank, skempi_support, score_composite)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, rows)
        self.conn.commit()

    def get_session_history(self, limit=20):
        cur = self.conn.execute("""
            SELECT session_id, input_name, input_target, created_at,
                   n_candidates, matched_pdb, runtime_sec, result_summary, status
            FROM analysis_sessions
            ORDER BY created_at DESC LIMIT ?
        """, (limit,))
        return [dict(r) for r in cur.fetchall()]

    def get_db_stats(self):
        return {
            'skempi_total':   self.conn.execute("SELECT COUNT(*) FROM skempi_mutations").fetchone()[0],
            'with_ddG':       self.conn.execute("SELECT COUNT(*) FROM skempi_mutations WHERE ddG_kcal IS NOT NULL").fetchone()[0],
            'with_koff':      self.conn.execute("SELECT COUNT(*) FROM skempi_mutations WHERE koff_ratio IS NOT NULL").fetchone()[0],
            'pdb_structures': self.conn.execute("SELECT COUNT(*) FROM pdb_structures").fetchone()[0],
            'residues':       self.conn.execute("SELECT COUNT(*) FROM pdb_residues").fetchone()[0],
            'interface_res':  self.conn.execute("SELECT COUNT(*) FROM pdb_residues WHERE is_interface=1").fetchone()[0],
            'total_sessions': self.conn.execute("SELECT COUNT(*) FROM analysis_sessions").fetchone()[0],
            'total_candidates': self.conn.execute("SELECT COUNT(*) FROM candidates").fetchone()[0],
        }
