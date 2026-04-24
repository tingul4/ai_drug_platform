"""
Build the PAI-1 mimicking peptide dataset (子計畫 3 Year-1 S1: Data Harmonization).

Source
------
Three FITC-Ahx labeled peptides synthesised by GenScript (Order U230SGAJG0),
COA PDFs in dataset/peptide_pai1/raw/.

Design (per document/agent/peptide_pipeline_plan.md §1 and §5):
* Naked FASTA for downstream computation (strip FITC-Ahx label).
* JSON metadata per peptide: batch_id, MW_theoretical, HPLC_purity, modification,
  source COA path — preserved so simulations stay traceable to the physical lot.
* MW sanity check: recomputed naked-peptide + FITC-Ahx adduct must match COA
  theoretical MW within 0.1 Da.

Outputs
-------
dataset/peptide_pai1/peptides.fasta   (naked sequences, no FITC-Ahx)
dataset/peptide_pai1/peptides.json    (per-peptide metadata + provenance)
dataset/peptide_pai1/schema.json      (JSON schema for peptides.json)
"""

import json
from pathlib import Path

ROOT = Path(__file__).parent.parent
OUT_DIR = ROOT / "dataset" / "peptide_pai1"
RAW_DIR = OUT_DIR / "raw"

# Monoisotopic-free, standard average amino-acid MW (Da).
AA_MW = {
    "A":  89.09, "R": 174.20, "N": 132.12, "D": 133.10, "C": 121.16,
    "E": 147.13, "Q": 146.15, "G":  75.07, "H": 155.16, "I": 131.17,
    "L": 131.17, "K": 146.19, "M": 149.21, "F": 165.19, "P": 115.13,
    "S": 105.09, "T": 119.12, "W": 204.22, "Y": 181.19, "V": 117.15,
}
WATER = 18.015
# FITC-Ahx N-terminal adduct (fluorescein-5-isothiocyanate + 6-aminohexanoic
# linker minus the water from amide bond formation). Calibrated against the
# three GenScript COA theoretical MWs (matches to 0.01 Da).
FITC_AHX_ADDUCT = 502.54

# Ground truth from the three GenScript Certificates of Analysis in RAW_DIR.
PEPTIDES = [
    {
        "id":              "pai1_peptide_69-80",
        "product_name":    "69-80",
        "sequence":        "KGMAPALRHLYK",
        "length":          12,
        "n_term_mod":      "FITC-Ahx",
        "mw_theoretical":  1887.24,
        "hplc_purity_pct": 98.0,
        "order_id":        "U230SGAJG0_9",
        "lot_no":          "U230SGAJG0-9/PE1018",
        "coa_path":        "raw/peptide_69-80_COA.pdf",
        "vendor":          "GenScript",
        "synthesis_date":  "2024-03-22",
    },
    {
        "id":              "pai1_peptide_76-88",
        "product_name":    "76-88",
        "sequence":        "RHLYKELMGPWNK",
        "length":          13,
        "n_term_mod":      "FITC-Ahx",
        "mw_theoretical":  2174.51,
        "hplc_purity_pct": 99.3,
        "order_id":        "U230SGAJG0_11",
        "lot_no":          "U230SGAJG0-11/PE1020",
        "coa_path":        "raw/peptide_76-88_COA.pdf",
        "vendor":          "GenScript",
        "synthesis_date":  "2024-03-22",
    },
    {
        "id":              "pai1_peptide_69-88",
        "product_name":    "69-88",
        "sequence":        "KGMAPALRHLYKELMGPWNK",
        "length":          20,
        "n_term_mod":      "FITC-Ahx",
        "mw_theoretical":  2843.36,
        "hplc_purity_pct": 98.2,
        "order_id":        "U230SGAJG0_7",
        "lot_no":          "U230SGAJG0-7/PE1016",
        "coa_path":        "raw/peptide_69-88_COA.pdf",
        "vendor":          "GenScript",
        "synthesis_date":  "2024-03-26",
    },
]

SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "title":   "PAI-1 mimicking peptide dataset (子計畫 3 S1)",
    "type":    "object",
    "required": ["target", "receptor_candidate", "peptides"],
    "properties": {
        "target":             {"type": "string"},
        "receptor_candidate": {"type": "string"},
        "source_document":    {"type": "string"},
        "peptides": {
            "type":  "array",
            "items": {
                "type": "object",
                "required": [
                    "id", "product_name", "sequence", "length", "n_term_mod",
                    "mw_theoretical", "mw_computed_naked", "mw_computed_with_mod",
                    "hplc_purity_pct", "order_id", "lot_no", "coa_path",
                ],
                "properties": {
                    "id":                   {"type": "string"},
                    "product_name":         {"type": "string"},
                    "sequence":             {"type": "string", "pattern": "^[ACDEFGHIKLMNPQRSTVWY]+$"},
                    "length":               {"type": "integer", "minimum": 1},
                    "n_term_mod":           {"type": ["string", "null"]},
                    "mw_theoretical":       {"type": "number"},
                    "mw_computed_naked":    {"type": "number"},
                    "mw_computed_with_mod": {"type": "number"},
                    "mw_delta":             {"type": "number"},
                    "hplc_purity_pct":      {"type": "number", "minimum": 0, "maximum": 100},
                    "order_id":             {"type": "string"},
                    "lot_no":               {"type": "string"},
                    "coa_path":             {"type": "string"},
                    "vendor":               {"type": "string"},
                    "synthesis_date":       {"type": "string"},
                },
            },
        },
    },
}


def naked_mw(seq: str) -> float:
    return sum(AA_MW[a] for a in seq) - (len(seq) - 1) * WATER


def build():
    if not RAW_DIR.exists():
        raise FileNotFoundError(f"Missing COA dir: {RAW_DIR}")

    records = []
    for p in PEPTIDES:
        seq = p["sequence"]
        if len(seq) != p["length"]:
            raise ValueError(f"{p['id']}: declared length {p['length']} != |seq| {len(seq)}")
        if (RAW_DIR / Path(p["coa_path"]).name).exists() is False:
            raise FileNotFoundError(f"{p['id']}: COA PDF not found: {p['coa_path']}")

        mw_naked = round(naked_mw(seq), 2)
        mw_with_mod = round(mw_naked + FITC_AHX_ADDUCT, 2) if p["n_term_mod"] == "FITC-Ahx" else mw_naked
        delta = round(mw_with_mod - p["mw_theoretical"], 3)
        if abs(delta) > 0.1:
            raise ValueError(
                f"{p['id']}: MW mismatch | computed {mw_with_mod:.2f} "
                f"vs COA {p['mw_theoretical']:.2f} (Δ={delta:+.3f})"
            )

        records.append({
            **p,
            "mw_computed_naked":    mw_naked,
            "mw_computed_with_mod": mw_with_mod,
            "mw_delta":             delta,
        })

    manifest = {
        "target":             "LRP1 CR cluster II (PDB 1J8E)",
        "receptor_candidate": "PDB:1J8E",
        "source_document":    "document/human/工程CM03-結合-unify-20260201-v35.pdf (p.1)",
        "note": (
            "PAI-1 mimicking peptides designed to compete with PAI-1 for LRP1 binding. "
            "FITC-Ahx is an experimental fluorescence label (for LRP1 binding assay) and "
            "is NOT part of the therapeutic moiety — computational modeling must use "
            "the naked sequence from 'sequence', not the labeled form."
        ),
        "peptides": records,
    }

    (OUT_DIR / "peptides.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (OUT_DIR / "schema.json").write_text(
        json.dumps(SCHEMA, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    fasta_lines = []
    for r in records:
        fasta_lines.append(
            f">{r['id']} | {r['product_name']} | {r['length']}AA | "
            f"lot={r['lot_no']} | naked_MW={r['mw_computed_naked']:.2f}"
        )
        fasta_lines.append(r["sequence"])
    (OUT_DIR / "peptides.fasta").write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")

    print(f"[build_peptide_dataset] wrote {len(records)} peptides → {OUT_DIR}")
    for r in records:
        print(
            f"  {r['id']:<22} {r['sequence']:<22} "
            f"naked={r['mw_computed_naked']:>7.2f}  +FITC-Ahx={r['mw_computed_with_mod']:>7.2f}  "
            f"COA={r['mw_theoretical']:>7.2f}  Δ={r['mw_delta']:+.3f}"
        )


if __name__ == "__main__":
    build()
