"""Subprocess helper for MHCflurry immunogenicity scoring.

Reads a JSON list of amino-acid sequences from stdin; writes a JSON dict
{sequence: {risk, risk_label}} to stdout.

Must be invoked with the .venv interpreter where MHCflurry 2.x is installed:
    echo '["KGMAPALRHLYK"]' | /raid/danielchen/ai_drug_platform/.venv/bin/python \
        /raid/danielchen/ai_drug_platform/engine/_mhcflurry_score.py
"""
from __future__ import annotations

import json
import sys
import warnings

from mhcflurry import Class1PresentationPredictor

# Taiwan Han-Chinese HLA-I panel (mirrors mhcflurry_immuno.py)
_HLA = [
    "HLA-A*11:01", "HLA-A*24:02", "HLA-A*33:03", "HLA-A*02:07", "HLA-A*02:01",
    "HLA-B*40:01", "HLA-B*46:01", "HLA-B*58:01", "HLA-B*13:01", "HLA-B*15:02",
    "HLA-C*01:02", "HLA-C*07:02", "HLA-C*08:01", "HLA-C*03:04", "HLA-C*03:02",
]
_WINDOW = 9
_ALLELE_BATCH = 6  # MHCflurry caps per-call allele count at 6


def _score_one(predictor, sequence: str) -> dict:
    if len(sequence) < _WINDOW:
        return {"risk": 0.0, "risk_label": "Low"}

    windows = [sequence[i : i + _WINDOW] for i in range(len(sequence) - _WINDOW + 1)]
    best: dict[int, float] = {}

    for start in range(0, len(_HLA), _ALLELE_BATCH):
        batch = _HLA[start : start + _ALLELE_BATCH]
        df = predictor.predict(peptides=windows, alleles=batch, verbose=0)
        for idx, (_, row) in enumerate(df.iterrows()):
            pct = float(row["presentation_percentile"])
            if idx not in best or pct < best[idx]:
                best[idx] = pct

    n = len(windows)
    n_strong = sum(1 for p in best.values() if p < 2.0)
    n_very = sum(1 for p in best.values() if p < 0.5)
    risk = round(0.7 * n_strong / n + 0.3 * n_very / n, 3)

    if risk >= 0.4:
        label = "High"
    elif risk >= 0.15:
        label = "Med"
    else:
        label = "Low"
    return {"risk": risk, "risk_label": label}


def main() -> None:
    sequences: list[str] = json.load(sys.stdin)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        predictor = Class1PresentationPredictor.load()
    results = {seq: _score_one(predictor, seq) for seq in sequences}
    print(json.dumps(results))


if __name__ == "__main__":
    main()
