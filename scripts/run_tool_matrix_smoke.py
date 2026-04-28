"""Run the 3 peptide x 3 Tool-A x 2 Tool-B smoke matrix and export results.

Default mode is fast smoke:
- Calls run_pipeline for each combination.
- Reports whether rho is available and whether alanine rows are real/stub.

Optional --real-ala mode is expensive:
- Forces combo-specific alanine sweep via run_alanine_scan_real before run_pipeline.
- Useful when you need rho for every Tool-A x Tool-B pair, not just cached ones.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from engine.analyzer import (
    PipelineTimeout,
    ProteinOptimizer,
    ToolAFailure,
    ToolBFailure,
)

PEPTIDES = [
    "KGMAPALRHLYK",
    "RHLYKELMGPWNK",
    "KGMAPALRHLYKELMGPWNK",
]
TOOL_AS = ["boltz2", "colabfold", "chai1"]
TOOL_BS = ["mmgbsa", "iptm"]


def parse_combo(text: str) -> tuple[str, str, str]:
    parts = text.split("::")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            "--skip-combo must be formatted as <sequence>::<tool_a>::<tool_b>"
        )
    sequence, tool_a, tool_b = [p.strip() for p in parts]
    if not sequence or not tool_a or not tool_b:
        raise argparse.ArgumentTypeError("--skip-combo fields cannot be empty")
    return sequence, tool_a, tool_b


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--output",
        default="/tmp/smoke_18_results.json",
        help="Output JSON path",
    )
    p.add_argument(
        "--real-ala",
        action="store_true",
        help="Force combo-specific run_alanine_scan_real before run_pipeline (slow).",
    )
    p.add_argument(
        "--force-live",
        action="store_true",
        help="Remove combo-relevant caches before each run so from_cache is avoided.",
    )
    p.add_argument(
        "--skip-combo",
        action="append",
        type=parse_combo,
        default=[],
        help="Skip one combo, format: <sequence>::<tool_a>::<tool_b>. Repeatable.",
    )
    return p.parse_args()


def _remove_if_exists(path: Path) -> bool:
    if not path.exists():
        return False
    if path.is_file():
        path.unlink()
        return True
    return False


def purge_combo_caches(
    opt: ProteinOptimizer, sequence: str, tool_a: str, tool_b: str
) -> int:
    """Delete cache files that could make this combo return from_cache=True."""
    removed = 0
    paths = [
        opt._cache_path(sequence, tool_a=tool_a),
        opt._cache_path(sequence),
        opt._ala_cache_path(sequence, tool_a=tool_a, tool_b=tool_b),
        opt._ala_cache_path(sequence, tool_a=tool_a),
        opt._ala_cache_path(sequence),
    ]
    for p in paths:
        if _remove_if_exists(p):
            removed += 1
    return removed


def main() -> int:
    args = parse_args()
    opt = ProteinOptimizer()
    rows = []
    skip = set(args.skip_combo)

    print(
        "sequence\ttool_a\ttool_b\tok\tfrom_cache\trho\treliability\ttool_a_s\ttool_b_s\tpost_s\ttotal_s\tbinding_energy"
    )

    for seq in PEPTIDES:
        for tool_a in TOOL_AS:
            for tool_b in TOOL_BS:
                combo = (seq, tool_a, tool_b)
                if combo in skip:
                    rec = {
                        "sequence": seq,
                        "tool_a": tool_a,
                        "tool_b": tool_b,
                        "ok": True,
                        "skipped": True,
                        "skip_reason": "requested by user",
                    }
                    rows.append(rec)
                    print(
                        "\t".join(
                            [
                                seq,
                                tool_a,
                                tool_b,
                                "True",
                                "-",
                                "-",
                                "-",
                                "-",
                                "-",
                                "-",
                                "-",
                                "SKIPPED",
                            ]
                        )
                    )
                    continue

                t0 = time.perf_counter()
                rec = {
                    "sequence": seq,
                    "tool_a": tool_a,
                    "tool_b": tool_b,
                    "ok": False,
                    "error": None,
                    "skipped": False,
                }
                try:
                    if args.force_live:
                        rec["removed_cache_files"] = purge_combo_caches(
                            opt, seq, tool_a, tool_b
                        )
                    if args.real_ala:
                        opt.run_alanine_scan_real(seq, tool_a=tool_a, tool_b=tool_b)
                    r = opt.run_pipeline(seq, tool_a=tool_a, tool_b=tool_b)
                    timings = r.get("timings", {})
                    rec.update(
                        {
                            "ok": True,
                            "from_cache": bool(r.get("from_cache")),
                            "rho": r.get("rho"),
                            "reliability": r.get("reliability"),
                            "binding_energy": r.get("binding_energy"),
                            "tool_a_s": timings.get("tool_a_s"),
                            "tool_b_s": timings.get("tool_b_s"),
                            "post_s": timings.get("post_s"),
                            "total_s": timings.get("total_s"),
                            "candidates_len": len(r.get("candidates") or []),
                            "alanine_len": len(r.get("alanine_scan") or []),
                            "alanine_has_real": any(
                                (
                                    x.get("predicted_ddg") is not None
                                    and x.get("mutant_dG") is not None
                                )
                                for x in (r.get("alanine_scan") or [])
                            ),
                        }
                    )
                except (ToolAFailure, ToolBFailure, PipelineTimeout, Exception) as e:
                    rec["error"] = f"{type(e).__name__}: {e}"

                rec["wall_s"] = round(time.perf_counter() - t0, 3)
                rows.append(rec)
                print(
                    "\t".join(
                        [
                            rec["sequence"],
                            rec["tool_a"],
                            rec["tool_b"],
                            str(rec["ok"]),
                            str(rec.get("from_cache")),
                            str(rec.get("rho")),
                            str(rec.get("reliability")),
                            str(rec.get("tool_a_s")),
                            str(rec.get("tool_b_s")),
                            str(rec.get("post_s")),
                            str(rec.get("total_s")),
                            str(rec.get("binding_energy")),
                        ]
                    )
                )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows, indent=2, ensure_ascii=False))
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
