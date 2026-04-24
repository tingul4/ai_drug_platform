"""
Mechanism-consistency benchmark: K→A alanine substitutions on PAI-1 seeds.

Reviewer-facing claim (張老師's Mechanistic Consistency):
    Stefansson et al. (full-length PAI-1 K→A mutants) show measurable loss of
    LRP1 affinity when K69 / K80 / K88 are mutated to Ala. A platform that can
    identify binding hotspots correctly must reproduce the *direction* of that
    trend on the 3 PAI-1 mimicking seeds — WT scores strongly, K→A scores
    weaker — regardless of the absolute metric. (Year-1 proxy is SKEMPI ΔΔG,
    not AF-Multimer ipTM; ipTM arrives with roadmap S3.)

Per-seed K coverage (by PAI-1 residue → seed position):
    69-80 seed (KGMAPALRHLYK):     K69→1, K80→12
    76-88 seed (RHLYKELMGPWNK):    K80→5, K88→13
    69-88 seed (KGMAPALRHLYK...):  K69→1, K80→12, K88→20

For each seed we report: WT, each K→A single, and the combined K→A multi-mutant
(2 or 3 mutations). All K residues are forced is_interface=True because
Stefansson identifies them as interface residues — this reflects the biological
claim rather than the position-based heuristic used by default alanine_scan.

Outputs:
    dataset/peptide_pai1/benchmark/ala_consistency.json
    dataset/peptide_pai1/benchmark/ala_consistency.md
"""

import json
import math
from pathlib import Path

from engine.analyzer import ProteinOptimizer

ROOT = Path(__file__).parent.parent
PEPTIDES_JSON = ROOT / "dataset" / "peptide_pai1" / "peptides.json"
OUT_DIR       = ROOT / "dataset" / "peptide_pai1" / "benchmark"

# PAI-1 absolute residue → (seed product_name, seed-local position)
# Only residues present in each seed are listed.
SEED_K_MAP = {
    "69-80": {"K69": 1,  "K80": 12},
    "76-88": {"K80": 5,  "K88": 13},
    "69-88": {"K69": 1,  "K80": 12, "K88": 20},
}


def _mutate(seq, pos_list, to='A'):
    """Apply 1-based point mutations, return (mutated_seq, mut_label_list)."""
    s = list(seq)
    labels = []
    for pos in pos_list:
        labels.append(f"{s[pos-1]}{pos}{to}")
        s[pos-1] = to
    return ''.join(s), labels


def _score_ala_mutations(model, seq, mutations):
    """
    Score a specific (wt_aa, pos, mut_aa) list under is_interface=True.
    Mirrors the aggregation rules in ProteinOptimizer.score_variants but
    skips the MHC / aggregation / solubility calls (affinity-only benchmark).
    """
    total_ddG = 0.0
    total_koff = 1.0
    total_support = 0
    per_site = []
    for (wt_aa, pos, mut_aa) in mutations:
        ddG, conf, n = model.predict_ddG(wt_aa, mut_aa, is_interface=True,
                                         is_buried=True, position=pos, sequence=seq)
        kr, kn = model.predict_koff(wt_aa, mut_aa)
        total_ddG += ddG
        total_koff *= kr
        total_support += n
        per_site.append({
            "pos":          pos,
            "wt_aa":        wt_aa,
            "mut_aa":       mut_aa,
            "ddG_kcal":     round(ddG, 3),
            "confidence":   round(conf, 3),
            "skempi_n":     n,
            "koff_ratio":   round(kr, 4),
        })

    # Kd fold change: ΔΔG = RT·ln(Kd_mut/Kd_WT), RT = 0.592 kcal/mol at 298 K
    # ⇒ Kd_mut / Kd_WT = 10^(ΔΔG / 1.363)
    pKd_shift = -total_ddG / 1.363
    kd_fold = math.pow(10.0, total_ddG / 1.363)

    return {
        "n_mutations":      len(mutations),
        "ddG_total":        round(total_ddG, 3),
        "pKd_shift":        round(pKd_shift, 3),
        "kd_fold_change":   round(kd_fold, 3),
        "koff_ratio_total": round(total_koff, 4),
        "skempi_n_total":   total_support,
        "per_site":         per_site,
    }


def build_benchmark_for_seed(opt, seed):
    name = seed['product_name']
    seq  = seed['sequence']
    k_map = SEED_K_MAP[name]

    rows = []
    # WT baseline
    rows.append({
        "label":      "WT",
        "mutations":  [],
        "sequence":   seq,
        **_score_ala_mutations(opt.model, seq, []),
    })

    # Singles
    for k_label, pos in k_map.items():
        mut_seq, mut_labels = _mutate(seq, [pos], to='A')
        scored = _score_ala_mutations(opt.model, seq, [(seq[pos-1], pos, 'A')])
        rows.append({
            "label":      f"{k_label}A",
            "mutations":  mut_labels,
            "sequence":   mut_seq,
            **scored,
        })

    # Multi-mutant: all Ks in this seed → A simultaneously
    if len(k_map) >= 2:
        positions = list(k_map.values())
        mut_seq, mut_labels = _mutate(seq, positions, to='A')
        triple_spec = [(seq[p-1], p, 'A') for p in positions]
        scored = _score_ala_mutations(opt.model, seq, triple_spec)
        rows.append({
            "label":      "+".join(f"{k}A" for k in k_map) + " (combined)",
            "mutations":  mut_labels,
            "sequence":   mut_seq,
            **scored,
        })

    return {
        "seed_id":      seed['id'],
        "product_name": name,
        "length":       seed['length'],
        "pai1_range":   name,  # seed product_name encodes PAI-1 residue range
        "wt_sequence":  seq,
        "k_sites":      k_map,
        "rows":         rows,
    }


def render_markdown(benchmark):
    lines = [
        "# Alanine Mechanism-Consistency Benchmark",
        "",
        "Year-1 baseline. Metric = SKEMPI-calibrated ΔΔG (not AF-Multimer ipTM — S3 pending).",
        "All K residues scored with `is_interface=True` per Stefansson 2004 PAI-1–LRP1 interface assignment.",
        "",
        "**Expected mechanism**: WT binding strongest; K→A mutations at interface hotspots weaken binding "
        "(ΔΔG > 0, Kd fold > 1, Δ pKd < 0, koff ratio > 1).",
        "",
    ]
    for seed in benchmark['seeds']:
        lines += [
            f"## Seed {seed['product_name']} — PAI-1 residues {seed['pai1_range']}",
            "",
            f"WT: `{seed['wt_sequence']}` ({seed['length']} AA). K sites in this seed: "
            + ", ".join(f"{k}→pos{v}" for k, v in seed['k_sites'].items()),
            "",
            "| Variant | Sequence | ΔΔG (kcal/mol) | Δ pKd | Kd_mut/Kd_WT | koff ratio | SKEMPI n |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
        for r in seed['rows']:
            ddg = r['ddG_total']
            lines.append(
                f"| {r['label']} | `{r['sequence']}` | "
                f"{ddg:+.3f} | {r['pKd_shift']:+.3f} | "
                f"{r['kd_fold_change']:.2f}× | {r['koff_ratio_total']:.3f} | {r['skempi_n_total']} |"
            )
        lines.append("")

    lines += [
        "## Interpretation",
        "",
        "- Positive ΔΔG and Kd fold > 1 on every K→A row = platform reproduces the Stefansson 2004 "
        "direction on all three PAI-1 mimicking peptides.",
        "- The `triple / combined` row on 69-88 is a superposition prediction: SKEMPIModel is trained "
        "on single-point ΔΔG, so multi-site sums are additive extrapolations (flagged in "
        "`analyzer.score_variants` when |ΔΔG| > 10 kcal/mol).",
        "",
        "## External-reference placeholder",
        "",
        "Fill when Stefansson paper Table values are extracted:",
        "",
        "| Mutant | Reported fold change (SPR K_D) | Our Kd_mut/Kd_WT | Direction match |",
        "|---|---|---|---|",
        "| K69A full-length PAI-1 | TBD | see above | — |",
        "| K80A full-length PAI-1 | TBD | see above | — |",
        "| K88A full-length PAI-1 | TBD | see above | — |",
    ]
    return "\n".join(lines) + "\n"


def build():
    manifest = json.loads(PEPTIDES_JSON.read_text())
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    opt = ProteinOptimizer()

    benchmark = {
        "metric":         "SKEMPI-calibrated ΔΔG (Year-1 proxy; ipTM pending roadmap S3)",
        "receptor":       manifest.get("receptor_candidate", "PDB:1J8E"),
        "reference":      "Stefansson et al. 2004, J Biol Chem — PAI-1 K→A mutants lose LRP1 affinity",
        "is_interface":   "forced True for all K sites (matches Stefansson interface assignment)",
        "seeds":          [],
    }
    for seed in manifest['peptides']:
        benchmark['seeds'].append(build_benchmark_for_seed(opt, seed))

    (OUT_DIR / "ala_consistency.json").write_text(
        json.dumps(benchmark, indent=2, ensure_ascii=False), encoding='utf-8'
    )
    (OUT_DIR / "ala_consistency.md").write_text(
        render_markdown(benchmark), encoding='utf-8'
    )

    # Console summary
    print(f"[benchmark] {len(benchmark['seeds'])} seeds → {OUT_DIR}")
    for seed in benchmark['seeds']:
        print(f"  {seed['product_name']}:")
        for r in seed['rows']:
            print(f"    {r['label']:<28}  ΔΔG={r['ddG_total']:+7.3f}  "
                  f"ΔpKd={r['pKd_shift']:+6.3f}  Kd_fold={r['kd_fold_change']:6.2f}×  "
                  f"koff={r['koff_ratio_total']:6.3f}")


if __name__ == "__main__":
    build()
