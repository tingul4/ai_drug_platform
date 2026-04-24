"""Small-molecule screening loader.

For the POC we expose the results that `dataset/kras_g12d_poc/` already produced
(Vina docking + ADMET-AI + koff proxy + PyMOO NSGA-II). The full screening
pipeline lives in the dataset folder as standalone scripts; the API layer here
just reads the TSV outputs and ships them to the frontend as JSON.

Extending later: wrap those scripts into `run_pipeline(target_pdb, smiles)` and
return its output in the same shape as `load_poc_results()`.
"""
from pathlib import Path
import pandas as pd

POC_DIR = Path(__file__).parent.parent / "dataset" / "kras_g12d_poc"


def _read_tsv(p: Path) -> list[dict]:
    if not p.exists():
        return []
    df = pd.read_csv(p, sep="\t")
    # pandas NaN -> None for JSON
    return [{k: (None if pd.isna(v) else v) for k, v in row.items()}
            for row in df.to_dict(orient="records")]


def load_poc_results() -> dict:
    """Return POC screening results in a frontend-friendly shape."""
    summary = _read_tsv(POC_DIR / "results" / "summary_table.tsv")
    ranked  = _read_tsv(POC_DIR / "results" / "candidates_ranked.tsv")
    front   = _read_tsv(POC_DIR / "results" / "pareto_front.tsv")

    # Build id -> pareto_rank map from ranked
    rank_by_id = {r["id"]: int(r.get("pareto_rank", 999)) for r in ranked}
    for c in summary:
        c["pareto_rank"] = rank_by_id.get(c["id"], 999)

    return {
        "target":      {"pdb_id": "7RPZ",
                         "uniprot": "P01116",
                         "name": "KRAS G12D",
                         "disease": "Pancreatic ductal adenocarcinoma (PDAC)"},
        "pipeline":    ["RDKit ETKDGv3 3D embed",
                         "Meeko PDBQT prep",
                         "AutoDock Vina 1.2.3 docking",
                         "ADMET-AI (Graph Attention Network)",
                         "Descriptor-based koff proxy",
                         "PyMOO NSGA-II / Non-dominated sorting"],
        "n_candidates":    len(summary),
        "n_pareto_front":  len(front),
        "candidates":      summary,
        "pareto_front_ids": [c["id"] for c in front],
    }


def poc_plot_path() -> Path:
    return POC_DIR / "results" / "pareto_plot.png"
