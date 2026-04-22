"""
Build a 9-mer MHC-II immunogenicity PSSM from the IEDB 2010 benchmark
(44,541 binding affinities, 28 alleles — http://tools.iedb.org/mhcii/download/).

Strategy
--------
* Binder   : IC50 <= 500 nM   (standard MHC-II binder threshold)
* Non-binder: IC50 >  500 nM
* For each labelled peptide we slide a 9-mer window over it and add each
  sub-9-mer to the binder / non-binder pool (weighted by 1 / num_windows
  per peptide so long peptides don't dominate).
* PSSM[pos][aa] = log2((p_binder + eps) / (p_background + eps))
* At inference a candidate sequence is scored by sliding 9-mer windows and
  taking the max window score, then normalised to 0-1 via binder-percentile
  anchors (5th = low risk, 95th = high risk).

Outputs:
  dataset/iedb_mhcii_2010/pssm.json   <- used by engine/analyzer.py
"""

import json, math, glob
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).parent.parent
DATA_DIR = ROOT / "dataset" / "iedb_mhcii_2010" / "class_II_all_split_5cv"
OUT_PATH = ROOT / "dataset" / "iedb_mhcii_2010" / "pssm.json"

AA = list("ACDEFGHIKLMNPQRSTVWY")
AA_IDX = {a: i for i, a in enumerate(AA)}
WINDOW = 9
BINDER_CUTOFF = 500.0
EPS = 1e-3

# Only use human HLA-DR/DP/DQ alleles (drop mouse) for human-immunogenicity PSSM
def is_human_allele(allele: str) -> bool:
    return allele.startswith("HLA-")

def load_records():
    """Yield (allele, peptide, ic50, is_binder)."""
    # Use only *_train_* files to avoid peeking at the canonical test split.
    files = sorted(glob.glob(str(DATA_DIR / "*_train_random_*.txt")))
    for fp in files:
        with open(fp) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 7:
                    continue
                species, allele, length, part, pep, ineq, ic50 = parts[:7]
                if not is_human_allele(allele):
                    continue
                try:
                    ic50v = float(ic50)
                except ValueError:
                    continue
                # Skip peptides containing non-standard residues
                pep = pep.upper().strip()
                if any(a not in AA_IDX for a in pep) or len(pep) < WINDOW:
                    continue
                yield allele, pep, ic50v, ic50v <= BINDER_CUTOFF

def build_pssm():
    binder_counts = [[0.0] * 20 for _ in range(WINDOW)]
    bg_counts     = [[0.0] * 20 for _ in range(WINDOW)]
    n_binder_pep = 0
    n_bg_pep     = 0
    binder_scores_raw = []   # used for percentile anchors, filled after PSSM built

    # First pass — accumulate position counts from all 9-mer sub-windows
    peptides = []
    for allele, pep, ic50, is_bind in load_records():
        peptides.append((pep, is_bind))
        n_wins = len(pep) - WINDOW + 1
        w = 1.0 / n_wins
        target = binder_counts if is_bind else bg_counts
        if is_bind:
            n_binder_pep += 1
        else:
            n_bg_pep += 1
        for i in range(n_wins):
            win = pep[i:i+WINDOW]
            for pos, aa in enumerate(win):
                target[pos][AA_IDX[aa]] += w

    print(f"[IEDB PSSM] human binder peptides = {n_binder_pep}, "
          f"non-binder peptides = {n_bg_pep}")

    # Normalise columns to probabilities, then log-ratio PSSM
    def col_probs(counts_row):
        tot = sum(counts_row) + 20 * EPS
        return [(c + EPS) / tot for c in counts_row]

    pssm = []
    for pos in range(WINDOW):
        pb = col_probs(binder_counts[pos])
        pg = col_probs(bg_counts[pos])
        pssm.append([math.log2(pb[i] / pg[i]) for i in range(20)])

    # Second pass — score every binder peptide to set percentile anchors
    def score_seq(seq):
        best = -1e9
        for i in range(len(seq) - WINDOW + 1):
            s = sum(pssm[pos][AA_IDX[seq[i+pos]]] for pos in range(WINDOW))
            if s > best:
                best = s
        return best

    binder_scores = [score_seq(p) for p, b in peptides if b]
    nonbinder_scores = [score_seq(p) for p, b in peptides if not b]
    binder_scores.sort()
    nonbinder_scores.sort()

    def pct(arr, q):
        if not arr: return 0.0
        k = max(0, min(len(arr)-1, int(q * (len(arr)-1))))
        return arr[k]

    anchors = {
        "binder_p05":     pct(binder_scores, 0.05),
        "binder_p50":     pct(binder_scores, 0.50),
        "binder_p95":     pct(binder_scores, 0.95),
        "nonbinder_p50":  pct(nonbinder_scores, 0.50),
        "nonbinder_p95":  pct(nonbinder_scores, 0.95),
    }

    out = {
        "source":       "IEDB MHC-II 2010 (classII_binding_data_Nov_16_2009)",
        "binder_cutoff_nM": BINDER_CUTOFF,
        "window":       WINDOW,
        "amino_acids":  AA,
        "n_binder_peptides": n_binder_pep,
        "n_nonbinder_peptides": n_bg_pep,
        "pssm":         pssm,       # 9 x 20 log2 ratios
        "anchors":      anchors,
    }
    OUT_PATH.write_text(json.dumps(out, indent=2))
    print(f"[IEDB PSSM] saved -> {OUT_PATH}")
    print(f"[IEDB PSSM] anchors: {anchors}")

if __name__ == "__main__":
    build_pssm()
