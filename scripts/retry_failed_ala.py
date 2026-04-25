"""Retry failed alanine-scan rows in a precomputed cache file.

Loads <seq>__ala.json, finds rows with `error`, re-runs Tool A + Tool B for
each, and writes successful results back in place. Skips rows that succeed
on retry; leaves the error in place if it fails again.
"""
import json, math, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from engine.analyzer import ProteinOptimizer, ToolAFailure, ToolBFailure, PipelineTimeout

SEQ = sys.argv[1] if len(sys.argv) > 1 else "KGMAPALRHLYKELMGPWNK"

opt = ProteinOptimizer()
cache_path = opt._ala_cache_path(SEQ)
rows = json.loads(cache_path.read_text())

wt = opt.load_cached(SEQ)
wt_dg = float(wt["binding_energy"])

n_failed = sum(1 for r in rows if r.get("error"))
print(f"[retry] {SEQ}: {n_failed}/{len(rows)} rows have errors", flush=True)

for idx, row in enumerate(rows):
    if not row.get("error"):
        continue
    i = row["position"] - 1
    mutant = SEQ[:i] + "A" + SEQ[i + 1:]
    print(f"[retry] pos {row['position']} {row['original']}->A ({mutant})", flush=True)
    try:
        mut_pdb = opt.run_tool_a(mutant, tool="boltz2")
        mut_dg = float(opt.run_tool_b(mut_pdb, sequence=mutant, tool="mmgbsa"))
    except (ToolAFailure, ToolBFailure, PipelineTimeout) as e:
        print(f"[retry]   FAIL again: {str(e)[:120]}", flush=True)
        continue

    ddg = mut_dg - wt_dg
    kd_fold = math.exp(ddg / 0.593)
    rows[idx] = {
        "position": i + 1,
        "original": SEQ[i],
        "mutant": "A",
        "predicted_ddg": round(ddg, 3),
        "kd_fold": round(kd_fold, 3),
        "mutant_dG": round(mut_dg, 3),
    }
    cache_path.write_text(json.dumps(rows, indent=2))
    print(f"[retry]   OK  ddG={ddg:+.3f}  kd_fold={kd_fold:.3g}", flush=True)

n_failed_after = sum(1 for r in rows if r.get("error"))
print(f"[retry] done: {n_failed_after}/{len(rows)} rows still failing", flush=True)
