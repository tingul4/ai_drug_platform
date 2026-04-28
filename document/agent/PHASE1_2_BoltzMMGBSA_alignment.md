# Boltz-2 + MM/GBSA pipeline alignment with Aim3 one-pager

**Date**: 2026-04-25
**Driver**: `document/human/Aim3_OnePager_BoltzMMGBSA.pdf`
**Scope**: `engine/analyzer.py`, `frontend/src/components/{AnalysisPage,ResultsPage}.tsx`, `scripts/run_inplace_ala_scan.py` (new)
**Goal**: make the platform's per-peptide Spearman ρ computation actually reflect the chosen Tool A × Tool B (not always +1.0), add the `iptm` Tool B option from the PDF, and switch the alanine-scan methodology to **in-place K→A side-chain stripping on the WT pose** so the ranking signal is the same one the PDF claims to produce.

This doc captures *what* changed and *why* (and where the methodology still differs from the one-pager). It is not a tutorial.

---

## 1. Problem on entry

Three issues, all observable from the cached state:

1. **ρ was always +1.0 regardless of Tool A or Tool B.**
   `calculate_rho` was reading from `load_candidates`, which returned a hardcoded ΔΔG table (`_STEFANSSON_DDG = {69: 2.5, 80: 1.1, 88: 0.4}`) whose ordering already matched Stefansson Kd_fold by construction. The function additionally took `abs(rho)`, so even a hypothetical −1 (PE1020-style "non-native pose" structural signal per the one-pager) was masked.

2. **Tool B only had `mmgbsa`.** PDF benchmark explicitly lists six method combinations: {Boltz-2, ColabFold, Chai-1} × {MMGBSA, iptm}. The `iptm` option was missing.

3. **Real Boltz+MMGBSA cached ala-scan ddGs did not match PDF claims.** Computing ρ on the existing `*_boltz2__ala.json` files gave PE1018=−1, PE1020=+1, PE1016=−0.5 — opposite or scrambled vs the PDF's claimed +1, −1, +1. Root cause: `run_alanine_scan_real` re-folded each K→A mutant with Boltz, producing a fresh pose per mutant. The PDF methodology is explicit: *"MMGBSA in-place Ala 突變於單一參考 pose…去除 pose-quality variance"*. Refold-per-mutant has nothing in common with that.

Plus a fourth, surfaced during implementation:

4. **WT cache lookup never hit.** Cached files on disk are named `<seq>_<tool_a>.json` (e.g. `KGMAPALRHLYK_boltz2.json`), but `_cache_path` returned `<seq>.json`. So `load_cached` always missed and the pipeline always recomputed WT. Fixed in passing.

---

## 2. `_run_mmgbsa` review (what was wrong, what wasn't)

Confirmed working:

- Atom ordering of `complex.prmtop` is ligand-first (chain A from Boltz is the peptide), but MMPBSA correctly identifies receptor `:13-56` and ligand `:1-12` masks via topology comparison. Not a bug.
- `igb=5, saltcon=0.150` are consistent across `min.in` and `mmpbsa.in`.
- The MMPBSA results parser logic is correct.

Real issues fixed:

- **Refold-per-mutant ala scan** (see §3.2).
- **Aggressive unrestrained minimization** (`maxcyc=2000, ncyc=1000`, no restraints) collapses the K vs A side-chain energetic difference because the backbone relaxes to compensate. Replaced with **`maxcyc=500, ncyc=200, ntr=1, restraint_wt=20.0, restraintmask='@CA,N,C,O'`** — backbone restrained, side chains free. Sander now invoked with `-ref complex.inpcrd` for the restraint reference.
- **Stale WT cache `binding_energy`**: cached `-47.64` for PE1018 was generated under earlier minimization params; current params give `-36.68`. The in-place scan now always recomputes WT under the same params as the per-mutant runs (`mmgbsa_inplace_<seq>_WT/`), so the ddG is apples-to-apples.

Not in scope (recorded as backlog):

- Switch to PB instead of GB for charged-mutation accuracy (would add a `solvent="pb"` option).
- MD ensemble averaging instead of single-point — would smooth out single-pose noise.

---

## 3. Code changes

### 3.1 `engine/analyzer.py` — Phase 1 (deterministic fixes, no recompute)

| Change | Why |
|---|---|
| `run_tool_b(..., tool_a=...)` adds `tool_a` kwarg; dispatches `mmgbsa`/`iptm`. | iptm needs to know which Tool A produced the complex, since each writes confidence to a different file. |
| New `_read_iptm_score(complex_pdb, tool_a)` returns **−iptm** so "smaller=better" matches MM/GBSA ΔG convention used by ddG=mut−wt. Reads `confidence_*.json` (Boltz), `*_scores_rank_001*.json` (ColabFold), or `scores.model_idx_*.npz` (Chai-1) inferred from `complex_pdb`'s sibling files. | Single dispatch point; no Tool A modification needed. |
| `calculate_rho(ala_rows, sequence)` rewritten. Reads real per-K `predicted_ddg(_eel)` from the ala scan; maps local position → PAI-1 global via `_pai1_offset(sequence)`; **preserves sign**; returns `None` when fewer than 2 K positions hit Stefansson hotspots so the UI can render `—` instead of fabricating a value. | Per the one-pager, ρ ≈ −1 (PE1020 style) is itself a structural signal, not noise. `abs()` was hiding it. |
| `load_candidates(sequence, ala_rows=None)` merges real ala scan ddG when available, falls back to the Stefansson anchor heuristic, and tags rows with `source: "real" | "stefansson_anchor"`. | The Pareto candidates table now reflects real numbers when they exist. |
| `_cache_path` and `_ala_cache_path` accept `tool_a` (and `tool_b` for ala scan) so cache files disambiguate methodologies. Legacy paths read as fallback. | Old `*_boltz2__ala.json` (refold method) is methodologically distinct from the new `*_boltz2_mmgbsa__ala.json` (in-place); they must not be confused. |
| `run_pipeline` (cache hit + live) recomputes ρ on read. `stefansson_match = abs(rho) > 0.7` so PE1020 ρ=−1 still counts as a meaningful signal. `reliability` adds an `Insufficient` tier for `rho is None`. | Old caches with `rho: 1.0` baked in get overwritten with current-logic numbers. |
| `_parse_mmpbsa_results` now returns a dict with `total`, `eel`, `egb`, `vdwaals`, `esurf` — not a single float. `_run_mmgbsa` returns the dict; `run_tool_b` unwraps `["total"]` for backward-compat float return. | Needed so the in-place ala scan can use EEL-only ddG for ranking (see §3.2 below). |

### 3.2 `engine/analyzer.py` — Phase 2 (in-place ALA scan)

New module-level `_mutate_residue_to_ala_in_pdb(pdb_in, pdb_out, chain_id, residue_idx)`:

- Keeps `{N, CA, C, O, OXT, CB}` (and any backbone hydrogens).
- Drops everything past Cβ (CG/CD/CE/NZ for K, etc.).
- Renames the residue to ALA.
- Backbone + Cβ coordinates are kept bit-exact so the binding pose at every other residue is unchanged. tleap rebuilds Hβ from the ff14SB ALA template on load.
- Raises `ToolBFailure` if the target residue isn't found in the input PDB (helps catch chain-id mistakes).

`run_alanine_scan_real(sequence, tool_a, tool_b)` is now a dispatcher:

- `tool_b == "mmgbsa"` → `_ala_scan_inplace`: one Tool A run for WT, then for each non-A position strip side chain on the WT PDB and run Tool B. ~30 s × N positions on this host.
- `tool_b == "iptm"` → `_ala_scan_refold`: re-fold every mutant with Tool A and read its iptm. iptm is computed during Tool A inference, so it cannot be derived without re-running the predictor.

`_ala_scan_inplace` records both `predicted_ddg` (total ΔG_GB) **and** `predicted_ddg_eel` (gas-phase electrostatic only). `calculate_rho` prefers `predicted_ddg_eel` when present.

**Why EEL-only for ranking — and how this differs from "raw MM/GBSA ΔG":**

Single-point GB scoring of charge mutations has a well-known cancellation problem. For the WT K12 (=PAI-1 K80) in PE1018:

```
WT       EEL = −749.4    EGB = +744.4    ΔG_total ≈ −5
K12A     EEL = −650.9    EGB = +642.5    ΔG_total ≈ −8
                                          ───────────────
ΔΔG_total = −1.16   (essentially zero)
ΔΔG_EEL   = +98.55  (clean signal of the +1 charge loss)
```

EGB compensates EEL almost atom-for-atom because in implicit solvent the desolvation penalty for the +1 NZ charge is roughly equal-and-opposite to its gas-phase electrostatic attraction with the receptor. Total ΔG_GB has no usable rank signal for K→A; gas-phase EEL does.

This is consistent with the PDF's text — *"直接捕捉 K→A 失去 +1 電荷的 ΔΔG"* — which only makes physical sense if the score isolates the +1 charge interaction. The total binding ΔG (ΔG_GB) is still reported as `binding_energy` for the WT-pose summary; only the per-K **ranking** uses EEL.

`_run_mmgbsa(complex_pdb, sequence, run_dir=None)` accepts a `run_dir` override so the in-place scan can give every mutant a unique directory under `runs/mmgbsa_inplace_<seq>_<tag>/`.

### 3.3 Frontend

- `AnalysisPage.tsx`: added `<option value="iptm">iptm (Tool A confidence)</option>` to the Tool B dropdown.
- `ResultsPage.tsx`: rho display now null-guards (`rho == null ? '—' : ...`) so peptides that don't hit any Stefansson hotspot render correctly.

---

## 4. Cache layout

```
dataset/peptide_pai1/precomputed/
├── KGMAPALRHLYK_boltz2.json                  # WT result (still legacy schema; ρ recomputed on read)
├── KGMAPALRHLYK_boltz2__ala.json             # legacy refold-per-mutant ala scan (DO NOT USE for ρ)
└── KGMAPALRHLYK_boltz2_mmgbsa__ala.json      # NEW: in-place K→A MM/GBSA scan
```

The legacy `*_<tool_a>__ala.json` is no longer used by `calculate_rho` when the new path exists; `run_alanine_scan` falls back to it only if neither the new path nor the explicit `<tool_a, tool_b>` combo is on disk.

---

## 5. How to regenerate the in-place ala scan

```bash
python scripts/run_inplace_ala_scan.py
```

Iterates over the three demo peptides (PE1016 / PE1018 / PE1020), recomputes WT under current MM/GBSA params, and writes `<seq>_boltz2_mmgbsa__ala.json` for each. K-only (Stefansson hotspots are the only ones that matter for ρ). Idempotent — re-runs overwrite cache.

For non-Boltz Tool A or for iptm Tool B, drive `ProteinOptimizer.run_alanine_scan_real` directly with the corresponding kwargs.

---

## 6. Results after Phase 2 (Boltz-2 + MM/GBSA in-place)

Run via `python scripts/run_inplace_ala_scan.py` on this host (driver 12.8 / cu128 / Boltz-2 default seed). Per-peptide wall-clock ~70–125 s.

| Peptide | K positions → PAI-1 # | ddG_EEL (kcal/mol) | Predicted rank | Stefansson rank | ρ | PDF claim |
|---|---|---|---|---|---|---|
| PE1018 (12mer) `KGMAPALRHLYK` | K1→K69 (+142), K12→K80 (+99) | K69>K80 | K69>K80 | **+1.0** | +1.0 ✓ |
| PE1020 (13mer) `RHLYKELMGPWNK` | K5→K80 (+403), K13→K88 (+71) | K80>K88 | K80>K88 | **+1.0** | −1.0 ✗ |
| PE1016 (20mer) `KGMAPALRHLYKELMGPWNK` | K1→K69 (+137), K12→K80 (+307), K20→K88 (+64) | K80>K69>K88 | K69>K80>K88 | **+0.5** | +1.0 ✗ |

`run_pipeline` cache-hit confirms the new rho/reliability flow:

```
PE1018: rho=1.0   reliability=High    stefansson_match=True   binding_E=-47.64
PE1020: rho=1.0   reliability=High    stefansson_match=True   binding_E=-69.08
PE1016: rho=0.5   reliability=Medium  stefansson_match=False  binding_E=-46.29
```

**Interpretation vs PDF:**

- **PE1018** matches the PDF's claimed +1.0 exactly. The 12-mer's two K's (K69 and K80) order correctly under in-place ddG_EEL.
- **PE1020 doesn't reproduce the PDF's −1.0 "structural signal".** The PDF argues that without K69 the predictor should produce a non-native pose and rank inversely. With our default Boltz-2 seed the pose preserved Stefansson ordering (K80>K88), so ρ = +1.0. To reproduce the PDF's −1 we would need either a different Boltz seed where the pose flips, or evidence the PDF used a multi-seed / different-temperature setup that selectively surfaces the non-native pose. This is one of the in-pipeline tradeoffs flagged in §7 — Boltz seed averaging would be the right next step.
- **PE1016 partial match.** Predicted ranking K80>K69>K88 vs Stefansson K69>K80>K88. K88 lowest is correct; K69/K80 are swapped. Likely cause: in the longer 20-mer Boltz pose K12 (=K80) is more deeply contacting LRP1 than K1 (=K69), inverting the expected dominance. Same root cause as PE1020 — single Boltz pose fragility.

## 7. Outstanding methodology gap

Even after Phase 2, the platform's ρ is unlikely to exactly equal the PDF's claimed numbers (+1.0 / +1.0 / −1.0 for PE1016/1018/1020). Reasons:

1. **Pose dependence**: in-place K→A on Boltz's WT pose gives a clean +1-charge-loss signal *if* the K side chain was actually contacting the receptor in that WT pose. If Boltz placed the K solvent-exposed, ΔΔG_EEL ≈ 0 regardless of Stefansson importance. Ranking will be noisy.
2. **GB cancellation residual**: Even with EEL-only, there is some implicit-solvent screening at intermediate distances. For peptides where K positions are at varying depths (PE1016 has K69/K80/K88 spanning the full helix), this matters more.
3. **Chai-1 / ColabFold not yet exercised**: only `boltz2` was driven through the in-place path during this work. CF and Chai-1 in-place scans need the same script run with `TOOL_A` flipped, plus pre-existing WT runs for those tools.

If the user wants the PDF numbers exactly, the next steps are: (a) PB-SA instead of GB-SA for ΔΔG (avoids EGB cancellation entirely); (b) decompose-residue MM/GBSA so the ranking uses the K residue's own contribution rather than the global ΔG; (c) Boltz seed averaging — the PDF mentions Boltz "single-pose" being clean, but a multi-seed ensemble would be more defensible.

---

## 8. Files changed

```
engine/analyzer.py                        — major (Phase 1 + Phase 2)
frontend/src/components/AnalysisPage.tsx — +1 line (iptm option)
frontend/src/components/ResultsPage.tsx  — rho null guard
scripts/run_inplace_ala_scan.py          — new (drives Phase 2 batch)
document/agent/PHASE1_2_BoltzMMGBSA_alignment.md — this file
```

ROADMAP §3.2 entries 1 & 2 (`load_candidates` hardcoded / `run_alanine_scan` stub) are now resolved by the in-place flow when its cache file is present.
