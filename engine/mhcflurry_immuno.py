"""
MHC-I immunogenicity risk model (子計畫 3 Year-1 S6, primary).

Wraps MHCflurry 2.0 `Class1PresentationPredictor` and aggregates 9-mer window
percentile-ranks across a Taiwanese Han-Chinese HLA-I panel. The panel choice
is the concrete implementation of the plan's 情境化校準 (Contextual Calibration)
requirement (工程CM03-v35 p.8): use Taiwan-population HLA frequencies so the
risk score reflects the patients actually enrolled at NCKU.

Aggregation rule per input peptide:
    for each 9-mer window:
        take min(percentile_rank) across the 15-allele panel     ← worst case
    risk_density  = fraction of windows where min_rank < 2.0     ← strong binders
    risk_density2 = fraction of windows where min_rank < 0.5     ← very strong
    risk          = 0.7 * risk_density + 0.3 * risk_density2 ∈ [0, 1]

References
----------
* O'Donnell et al., Cell Systems 2020 — MHCflurry 2.0
* Allele Frequency Net Database — Taiwan Han-Chinese HLA-A/B/C frequencies
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import List, Optional

# Taiwan Han-Chinese high-frequency HLA-I panel (allele, population frequency).
# Source: Allele Frequency Net DB + TWB (Taiwan Biobank) public summaries.
# Coverage: ~95% of Taiwan population covered by at least one allele in the panel.
HLA_I_TAIWAN_PANEL: List[str] = [
    # HLA-A
    "HLA-A*11:01",  # ~0.37
    "HLA-A*24:02",  # ~0.18
    "HLA-A*33:03",  # ~0.09
    "HLA-A*02:07",  # ~0.07
    "HLA-A*02:01",  # ~0.05
    # HLA-B
    "HLA-B*40:01",  # ~0.17
    "HLA-B*46:01",  # ~0.12
    "HLA-B*58:01",  # ~0.07
    "HLA-B*13:01",  # ~0.06
    "HLA-B*15:02",  # ~0.04
    # HLA-C
    "HLA-C*01:02",  # ~0.21
    "HLA-C*07:02",  # ~0.14
    "HLA-C*08:01",  # ~0.12
    "HLA-C*03:04",  # ~0.08
    "HLA-C*03:02",  # ~0.05
]

MHC_I_WINDOW = 9
STRONG_BINDER_PCT = 2.0
VERY_STRONG_PCT = 0.5


@dataclass
class MHCIRiskResult:
    risk: float                 # 0-1
    risk_density: float         # fraction windows with any allele < 2%
    risk_density2: float        # fraction windows with any allele < 0.5%
    n_windows: int
    top_hits: list              # up to 5 strongest (position, window, allele, percentile)
    strong_peptides: frozenset  # all 9-mer sequences with %rank < 2 (for Δ-from-WT Penalty)
    min_percentile: Optional[float]
    mean_percentile: Optional[float]


class MHCflurryImmunoModel:
    """MHCflurry-backed MHC-I immunogenicity scorer with a Taiwan HLA-I panel."""

    def __init__(self, alleles: Optional[List[str]] = None):
        from mhcflurry import Class1PresentationPredictor
        self.alleles = list(alleles) if alleles else list(HLA_I_TAIWAN_PANEL)
        # Suppress the chatty TF/torch device-detection warnings on load.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.predictor = Class1PresentationPredictor.load()
        self.source = f"MHCflurry 2.x, Taiwan HLA-I panel ({len(self.alleles)} alleles)"

    def score_sequence(self, sequence: str) -> Optional[MHCIRiskResult]:
        if len(sequence) < MHC_I_WINDOW:
            return None
        windows = [
            (i + 1, sequence[i:i + MHC_I_WINDOW])
            for i in range(len(sequence) - MHC_I_WINDOW + 1)
        ]
        peptides = [w[1] for w in windows]

        # MHCflurry caps `predict(alleles=[...])` at 6 alleles (one HLA genotype).
        # Our 15-allele Taiwan panel needs to be run in ≤6-sized batches; we keep
        # the per-peptide minimum presentation_percentile across all batches.
        BATCH = 6
        best_per_peptide: dict = {}
        for start in range(0, len(self.alleles), BATCH):
            batch = self.alleles[start:start + BATCH]
            df = self.predictor.predict(peptides=peptides, alleles=batch, verbose=0)
            for i, (_, row) in enumerate(df.iterrows()):
                pct = float(row["presentation_percentile"])
                prev = best_per_peptide.get(i)
                if prev is None or pct < prev["presentation_percentile"]:
                    best_per_peptide[i] = {
                        "best_allele":             str(row["best_allele"]),
                        "presentation_percentile": pct,
                        "affinity_nM":             float(row["affinity"]),
                    }

        per_window = []
        for i, (pos, win) in enumerate(windows):
            rec = best_per_peptide[i]
            per_window.append({
                "pos":                     pos,
                "peptide":                 win,
                "best_allele":             rec["best_allele"],
                "presentation_percentile": rec["presentation_percentile"],
                "affinity_nM":             rec["affinity_nM"],
            })

        n = len(per_window)
        hits_strong = [w for w in per_window if w["presentation_percentile"] < STRONG_BINDER_PCT]
        hits_very = [w for w in per_window if w["presentation_percentile"] < VERY_STRONG_PCT]
        risk_density = len(hits_strong) / n
        risk_density2 = len(hits_very) / n
        risk = round(0.7 * risk_density + 0.3 * risk_density2, 3)

        top_hits = sorted(per_window, key=lambda w: w["presentation_percentile"])[:5]
        min_pct = top_hits[0]["presentation_percentile"] if top_hits else None
        mean_pct = sum(w["presentation_percentile"] for w in per_window) / n

        return MHCIRiskResult(
            risk=risk,
            risk_density=round(risk_density, 3),
            risk_density2=round(risk_density2, 3),
            n_windows=n,
            top_hits=[
                (w["pos"], w["peptide"], w["best_allele"], round(w["presentation_percentile"], 3))
                for w in top_hits
            ],
            strong_peptides=frozenset(w["peptide"] for w in hits_strong),
            min_percentile=round(min_pct, 3) if min_pct is not None else None,
            mean_percentile=round(mean_pct, 3),
        )
