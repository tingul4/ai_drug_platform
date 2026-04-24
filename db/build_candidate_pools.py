"""
Build Year-1 candidate pools for the 3 PAI-1 mimicking peptide seeds (S8 runner).

Per document/agent/peptide_pipeline_plan.md + document/note.md alignment,
produces 4 pools total:
  - 3 independent pools: {69-80, 76-88, 69-88}, each ≥ 300 variants
  - 1 merged pool:       top-K by J(x) from each independent pool, re-ranked jointly

J(x) Year-1 proxies (see engine/analyzer.ProteinOptimizer.composite_Jx):
  f_bind  ← SKEMPI ΔΔG       (replaced by docking ΔG when S3/S4 land)
  f_admet ← heuristics       (replaced by ADMETlab 2.0 when S7 lands)
  f_synth ← 0.5 constant     (replaced by AiZynthFinder when S7 lands)
  Penalty ← MHC-I %rank<2 OR aggregation>0.4 (Tox21 flag arrives in S7)

Outputs:
  dataset/peptide_pai1/candidates/{69-80,76-88,69-88,merged}/variants.json
  dataset/peptide_pai1/candidates/summary.json
"""

import json
import time
from pathlib import Path

from engine.analyzer import ProteinOptimizer

ROOT = Path(__file__).parent.parent
PEPTIDES_JSON = ROOT / "dataset" / "peptide_pai1" / "peptides.json"
OUT_DIR       = ROOT / "dataset" / "peptide_pai1" / "candidates"

N_VARIANTS_PER_SEED    = 320    # roadmap S8: ≥ 300 per seed
MERGED_TOP_K_PER_SEED  = 100    # top-J(x) per seed fed into merged pool

# Keys that survive into JSON output (drops internal / non-serializable ones).
KEEP_KEYS = (
    'variant_id', 'mutations', 'sequence', 'is_wt', 'seed_origin',
    'ddG_binding', 'ddG_stability', 'pKd_shift', 'pKd_abs',
    'koff_relative', 'residence_time',
    'immunogenicity', 'mhc_i_introduces_new_binder', 'mhc_i_new_peptides',
    'aggregation', 'solubility', 'pI', 'mw_kDa', 'hbonds',
    'skempi_support', 'pareto_rank', 'warnings', 'Jx',
)


def _strip(v):
    return {k: v[k] for k in KEEP_KEYS if k in v}


def _pool_stats(scored):
    non_wt = [v for v in scored if not v.get('is_wt')]
    if not non_wt:
        return {}
    js = [v['Jx']['J'] for v in non_wt]
    pens = sum(1 for v in non_wt if v['Jx']['penalty'] > 0)
    n_p1 = sum(1 for v in non_wt if v.get('pareto_rank') == 1)
    return {
        'n_variants':      len(non_wt),
        'n_pareto_rank1':  n_p1,
        'n_penalty':       pens,
        'J_min':           round(min(js), 4),
        'J_max':           round(max(js), 4),
        'J_mean':          round(sum(js) / len(js), 4),
    }


def build_seed_pool(opt, seed, n_variants):
    seq = seed['sequence']
    t0 = time.time()
    scan = opt.alanine_scan(seq)
    variants = opt.generate_variants(seq, scan, max_variants=n_variants, strategy='mixed')
    scored = opt.score_variants(variants, scan)
    scored = opt.pareto_optimize(scored)
    for v in scored:
        v['seed_origin'] = seed['product_name']
        v['Jx'] = opt.composite_Jx(v)
    
    # Primary sort: Pareto rank (lower is better), secondary: J(x) (higher is better)
    wt = [v for v in scored if v.get('is_wt')]
    non_wt = sorted(
        [v for v in scored if not v.get('is_wt')],
        key=lambda v: (v['pareto_rank'], -v['Jx']['J'])
    )
    ordered = wt + non_wt
    return ordered, scan, round(time.time() - t0, 2)


def build_merged_pool(per_seed_pools, opt, top_k=MERGED_TOP_K_PER_SEED):
    """Take top-K (by J(x)) non-WT variants from each seed pool and re-rank jointly.
    
    Includes cross-seed de-duplication by sequence to ensure unique candidates.
    """
    combined_raw = []
    for pool in per_seed_pools:
        non_wt = [v for v in pool if not v.get('is_wt')]
        # Sort by Jx to pick top-K from this seed
        non_wt.sort(key=lambda v: -v['Jx']['J'])
        combined_raw.extend(non_wt[:top_k])

    # De-duplicate by sequence: if multiple seeds produce the same sequence, 
    # keep the one with the higher J(x).
    unique_variants = {}
    for v in combined_raw:
        seq = v['sequence']
        if seq not in unique_variants or v['Jx']['J'] > unique_variants[seq]['Jx']['J']:
            unique_variants[seq] = v
    
    combined = list(unique_variants.values())

    # Re-run Pareto on the merged set so ranks reflect the joint front.
    combined = opt.pareto_optimize(combined)
    # Final Sort: Pareto Rank 1st, J(x) 2nd
    combined.sort(key=lambda v: (v['pareto_rank'], -v['Jx']['J']))
    
    # Re-assign IDs for the merged pool to be clean
    for i, v in enumerate(combined):
        v['variant_id'] = f"M{i+1:04d}"
        
    return combined


def write_pool(pool_name, ordered, scan, seed_meta, runtime):
    pdir = OUT_DIR / pool_name
    pdir.mkdir(parents=True, exist_ok=True)
    payload = {
        'pool_name':      pool_name,
        'seed':           seed_meta,      # None for merged
        'runtime_sec':    runtime,
        'stats':          _pool_stats(ordered),
        'alanine_scan':   scan,           # None for merged
        'variants':       [_strip(v) for v in ordered],
    }
    (pdir / "variants.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=str),
        encoding='utf-8'
    )
    return payload['stats']


def build():
    manifest = json.loads(PEPTIDES_JSON.read_text())
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    opt = ProteinOptimizer()
    summary = {
        'generated_at': time.strftime('%Y-%m-%d %H:%M:%S'),
        'n_variants_per_seed': N_VARIANTS_PER_SEED,
        'merged_top_k_per_seed': MERGED_TOP_K_PER_SEED,
        'weights': {
            'w_bind':    opt.JX_W_BIND,
            'w_admet':   opt.JX_W_ADMET,
            'w_synth':   opt.JX_W_SYNTH,
            'w_penalty': opt.JX_W_PENALTY,
        },
        'proxy_sources': {
            'f_bind':  'SKEMPIModel ΔΔG (Year-1 proxy; S3/S4 → docking ΔG)',
            'f_admet': 'solubility/aggregation heuristic (Year-1 proxy; S7 → ADMETlab 2.0)',
            'f_synth': 'constant 0.5 (Year-1 placeholder; S7 → AiZynthFinder)',
            'penalty': 'new MHC-I 9-mer %rank<2 vs WT OR aggregation>0.4 (Tox21 pending S7)',
        },
        'pools': {},
    }

    per_seed_pools = []
    for seed in manifest['peptides']:
        name = seed['product_name']
        print(f"[pool] {name}: generating {N_VARIANTS_PER_SEED} variants...")
        ordered, scan, runtime = build_seed_pool(opt, seed, N_VARIANTS_PER_SEED)
        stats = write_pool(name, ordered, scan, seed, runtime)
        summary['pools'][name] = {**stats, 'runtime_sec': runtime}
        per_seed_pools.append(ordered)
        print(f"       {stats} in {runtime}s")

    print(f"[pool] merged: combining top-{MERGED_TOP_K_PER_SEED} from each seed...")
    t0 = time.time()
    merged = build_merged_pool(per_seed_pools, opt)
    # attach J(x) was already on each variant; Pareto updated above
    merged_runtime = round(time.time() - t0, 2)
    stats = write_pool('merged', merged, None, None, merged_runtime)
    summary['pools']['merged'] = {**stats, 'runtime_sec': merged_runtime}

    (OUT_DIR / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding='utf-8'
    )
    print(f"[pool] summary → {OUT_DIR / 'summary.json'}")
    print(json.dumps(summary['pools'], indent=2))


if __name__ == "__main__":
    build()
