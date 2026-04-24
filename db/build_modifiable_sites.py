"""
Build modifiable-site lists for the 3 PAI-1 mimicking peptides (roadmap S5).

Schema (per document/agent/roadmap.md S5):
    residue_index              — 1-based position in the seed sequence
    wt_aa                      — wild-type residue
    allowed_mutations          — ranked list of candidate substitutions
    interface_importance_score — ΔΔG_Ala from SKEMPIModel (kcal/mol; ≈ hotspot strength)
    interaction_type           — Year-1 proxy from AA class; S3/S4 replaces with PLIF fingerprint
    evidence_source            — SKEMPI support count + hotspot tag

Year-1 note: alanine_scan falls back to position-based interface heuristic when
no PDB is bound; the LRP1 CR7 receptor we prepared in S2 is not yet docked,
so `is_interface` is a sequence-position proxy. Upgrade path: S3 AF-Multimer
→ S4 PLIF to tag interaction_type with real donor/acceptor/salt-bridge flags.
"""

import json
from pathlib import Path

from engine.analyzer import AA_PROPS, ProteinOptimizer, blosum_score

ROOT = Path(__file__).parent.parent
PEPTIDES_JSON = ROOT / "dataset" / "peptide_pai1" / "peptides.json"
OUT_DIR = ROOT / "dataset" / "peptide_pai1" / "modifiable_sites"

# AA physicochemical class tags (Year-1 proxy for PLIF interaction_type).
AA_CLASS = {
    **{a: 'hydrophobic' for a in 'AVLIMP'},
    **{a: 'polar'       for a in 'STNQ'},
    **{a: 'charged_pos' for a in 'KRH'},
    **{a: 'charged_neg' for a in 'DE'},
    **{a: 'aromatic'    for a in 'FYW'},
    'C': 'thiol', 'G': 'flex',
}

MAX_ALLOWED_PER_SITE = 8


def _rank_allowed_mutations(opt, wt_aa):
    """Rank substitutions by SKEMPI improvement + BLOSUM conservativeness."""
    scored = []
    for aa in AA_PROPS:
        if aa == wt_aa:
            continue
        skempi_mean = opt.model.sub_mean.get((wt_aa, aa))       # lower = stronger binding
        skempi_n    = len(opt.model.sub_stats.get((wt_aa, aa), []))
        blosum      = blosum_score(wt_aa, aa)
        # Lower composite = better candidate.
        composite = (skempi_mean if skempi_mean is not None else 0.0) - 0.1 * blosum
        scored.append({
            'aa':          aa,
            'aa_class':    AA_CLASS.get(aa, 'other'),
            'skempi_ddG':  round(skempi_mean, 3) if skempi_mean is not None else None,
            'skempi_n':    skempi_n,
            'blosum':      blosum,
            '_composite':  composite,
        })
    scored.sort(key=lambda r: r['_composite'])
    out = []
    for r in scored[:MAX_ALLOWED_PER_SITE]:
        r.pop('_composite')
        out.append(r)
    return out


def build_sites_for_seed(opt, seed):
    seq = seed['sequence']
    scan = opt.alanine_scan(seq)
    scan_by_pos = {r['position']: r for r in scan}

    sites = []
    for pos in range(1, len(seq) + 1):
        wt_aa = seq[pos - 1]
        r = scan_by_pos.get(pos)
        if r is None:
            # Alanine residues are skipped by alanine_scan (wt == 'A'); report
            # them as unmodifiable-by-default so schema stays dense.
            sites.append({
                'residue_index':              pos,
                'wt_aa':                      wt_aa,
                'allowed_mutations':          [],
                'interface_importance_score': 0.0,
                'interaction_type':           AA_CLASS.get(wt_aa, 'other'),
                'evidence_source':            'alanine_scan_skipped (wt=A)',
                'modifiable':                 False,
            })
            continue
        modifiable = r['contribution'] != 'neutral' or r['is_interface']
        sites.append({
            'residue_index':              pos,
            'wt_aa':                      wt_aa,
            'allowed_mutations':          _rank_allowed_mutations(opt, wt_aa) if modifiable else [],
            'interface_importance_score': r['ddG_ala'],
            'interaction_type':           AA_CLASS.get(wt_aa, 'other'),
            'evidence_source': (
                f"alanine_scan: contribution={r['contribution']}, "
                f"is_interface={r['is_interface']}, skempi_n={r['skempi_n']}, "
                f"confidence={r['confidence']}"
            ),
            'modifiable': modifiable,
        })
    return sites


SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "title":   "Modifiable sites list (子計畫 3 S5)",
    "type":    "object",
    "required": ["seed_id", "sequence", "sites"],
    "properties": {
        "seed_id":  {"type": "string"},
        "sequence": {"type": "string", "pattern": "^[ACDEFGHIKLMNPQRSTVWY]+$"},
        "sites": {
            "type":  "array",
            "items": {
                "type": "object",
                "required": [
                    "residue_index", "wt_aa", "allowed_mutations",
                    "interface_importance_score", "interaction_type",
                    "evidence_source", "modifiable",
                ],
            },
        },
    },
}


def build():
    manifest = json.loads(PEPTIDES_JSON.read_text())
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "schema.json").write_text(
        json.dumps(SCHEMA, indent=2, ensure_ascii=False), encoding='utf-8'
    )

    opt = ProteinOptimizer()
    for seed in manifest['peptides']:
        sites = build_sites_for_seed(opt, seed)
        n_modifiable = sum(1 for s in sites if s['modifiable'])
        payload = {
            'seed_id':            seed['id'],
            'product_name':       seed['product_name'],
            'sequence':           seed['sequence'],
            'length':             seed['length'],
            'receptor_candidate': manifest['receptor_candidate'],
            'note':               (
                "Year-1 baseline. interaction_type is an AA-class proxy; S3 AF-Multimer "
                "+ S4 PLIF will replace it with real donor/acceptor/salt-bridge tags. "
                "is_interface inside evidence_source is a sequence-position heuristic "
                "until docking is wired."
            ),
            'n_sites':            len(sites),
            'n_modifiable':       n_modifiable,
            'sites':              sites,
        }
        out = OUT_DIR / f"{seed['product_name']}.json"
        out.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding='utf-8')
        print(f"[S5] {seed['product_name']:>6}: {n_modifiable}/{len(sites)} modifiable → {out.name}")


if __name__ == "__main__":
    build()
