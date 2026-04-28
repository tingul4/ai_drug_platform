"""Run the in-place K→A MM/GBSA scan for the three Aim3 demo peptides.

Drives `ProteinOptimizer._ala_scan_inplace` for each peptide: WT pose is
re-evaluated under current MM/GBSA params, then each non-A position is
side-chain-stripped on the WT pose and scored single-point. Output goes to
`dataset/peptide_pai1/precomputed/<seq>_<tool_a>_mmgbsa__ala.json`.

Per-peptide cost: ~30 s × N positions ≈ 5–10 min on this host. Safe to run
sequentially; safe to re-run (idempotent — overwrites cache file).
"""
import json
import shutil
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from engine.analyzer import (  # noqa: E402
    ProteinOptimizer,
    _mutate_residue_to_ala_in_pdb,
    PEPTIDE_CHAIN,
    RUNS_DIR,
    ToolBFailure,
    PipelineTimeout,
)
from scipy.stats import spearmanr  # noqa: E402

PEPTIDES = ["KGMAPALRHLYK", "RHLYKELMGPWNK", "KGMAPALRHLYKELMGPWNK"]
TOOL_A = "boltz2"
KD_BY_GLOBAL = ProteinOptimizer._STEFANSSON_KD_BY_GLOBAL


def scan_one(opt, seq):
    print(f"\n=== {seq} ===", flush=True)
    t0 = time.perf_counter()
    wt = opt.load_cached(seq, tool_a=TOOL_A)
    if wt is None:
        raise SystemExit(f"WT cache missing for {seq}; run the WT pipeline first.")
    wt_complex = wt["complex_pdb"]
    if not Path(wt_complex).exists():
        raise SystemExit(f"WT complex PDB does not exist on disk: {wt_complex}")

    # Recompute WT components under current params.
    wt_run = RUNS_DIR / f"mmgbsa_inplace_{seq}_WT"
    wt_run.mkdir(parents=True, exist_ok=True)
    shutil.copy(wt_complex, wt_run / "complex_mut.pdb")
    wt_components = opt._run_mmgbsa(str(wt_run / "complex_mut.pdb"), sequence=seq, run_dir=wt_run)
    wt_total = float(wt_components["total"])
    wt_eel = float(wt_components.get("eel", 0.0))
    print(f"WT total={wt_total:.3f}  EEL={wt_eel:.3f}", flush=True)

    rows = []
    positions = [(i, aa) for i, aa in enumerate(seq) if aa == "K"]  # K-only is enough for ρ
    for i, aa in positions:
        local = i + 1
        tag = f"K{local}A"
        run_dir = RUNS_DIR / f"mmgbsa_inplace_{seq}_{tag}"
        run_dir.mkdir(parents=True, exist_ok=True)
        mut_pdb = run_dir / "complex_mut.pdb"
        mutant_label = seq[:i] + "A" + seq[i + 1:]
        try:
            _mutate_residue_to_ala_in_pdb(
                wt_complex, str(mut_pdb),
                chain_id=PEPTIDE_CHAIN, residue_idx=local,
            )
            mut_components = opt._run_mmgbsa(str(mut_pdb), sequence=mutant_label, run_dir=run_dir)
        except (ToolBFailure, PipelineTimeout) as e:
            print(f"  K{local}A FAIL: {e}", flush=True)
            rows.append({"position": local, "original": "K", "mutant": "A",
                         "error": str(e)[:200], "method": "inplace"})
            continue
        mut_total = float(mut_components["total"])
        mut_eel = float(mut_components.get("eel", 0.0))
        ddg = mut_total - wt_total
        ddg_eel = mut_eel - wt_eel
        rows.append({
            "position":          local,
            "original":          "K",
            "mutant":            "A",
            "predicted_ddg":     round(ddg, 3),
            "predicted_ddg_eel": round(ddg_eel, 3),
            "kd_fold":           round(float(__import__("math").exp(ddg / 0.593)), 3),
            "mutant_dG":         round(mut_total, 3),
            "method":            "inplace",
        })
        print(f"  K{local}A  total={mut_total:.3f}  ddG={ddg:+.3f}  ddG_EEL={ddg_eel:+.3f}", flush=True)

    cache_path = opt._ala_cache_path(seq, tool_a=TOOL_A, tool_b="mmgbsa")
    cache_path.write_text(json.dumps(rows, indent=2))
    print(f"  -> wrote {cache_path.name}", flush=True)

    rho = opt.calculate_rho(rows, seq)
    elapsed = time.perf_counter() - t0
    print(f"  signed ρ vs Stefansson = {rho}  (elapsed {elapsed:.1f}s)", flush=True)
    return rho


def main():
    opt = ProteinOptimizer()
    rhos = {}
    for seq in PEPTIDES:
        rhos[seq] = scan_one(opt, seq)
    print("\n=== Summary ===")
    for seq, rho in rhos.items():
        print(f"  {seq}: ρ = {rho}")


if __name__ == "__main__":
    main()
