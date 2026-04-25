# AIM3 — PAI-1 Mimicking Peptide Optimization Platform

胜肽藥物優化平台：以 Boltz-2 / ColabFold / Chai-1（Tool A）做 peptide-receptor 複合物 3D 預測，再以 AmberTools MM/GBSA（Tool B, igb=5）算結合能。針對 PAI-1 mimicking peptides 對 LRP1 CR cluster II 的競爭結合進行設計與校準。

---

## Architecture

```
ai_drug_platform/
├── api/server.py              # Flask + 非同步 Job Queue
├── engine/analyzer.py         # 主 pipeline (Tool A → Tool B → ρ + Pareto)
├── db/build_lrp1_receptor.py  # 受體準備 (PDB 1J8E → 清洗 + 最小化)
├── dataset/
│   ├── lrp1_cr2/              # 已準備好的 LRP1 CR7 受體
│   └── peptide_pai1/
│       ├── peptides.json      # 三條 demo 胜肽元數據
│       └── precomputed/       # demo 胜肽的 cached pipeline 結果
├── third_party/               # 所有 conda env (見下節)
│   ├── gmx_env/               # 主 env: GROMACS + AmberTools + Flask
│   ├── boltz_env/             # Tool A: Boltz-2
│   ├── colabfold_env/         # Tool A: ColabFold
│   └── chai_env/              # Tool A: Chai-1
├── frontend/                  # Vite + React 前端
├── runs/                      # 每次 pipeline 的工作目錄 (gitignored)
└── requirements.txt           # gmx_env Python 套件 (其他 env 列在下節)
```

---

## Environment setup

平台採用 **multi-conda-env** 架構：四個工具因 `numpy` / `biopython` 等版本相互衝突，無法同 env 共存；每個工具一個獨立 conda env，analyzer 透過 subprocess 呼叫各 env 的 CLI。

### Prerequisites
- conda（建議 24+，需要 `--solver=libmamba`）
- NVIDIA GPU + driver ≥ 12.8（受體準備 + Boltz/ColabFold/Chai 推理用）

### 1. 主 env (`gmx_env`) — Flask 後端 + Tool B (MM/GBSA)

```bash
conda create -p ./third_party/gmx_env -c conda-forge -c bioconda \
    python=3.10 gromacs ambertools gmx_mmpbsa openmm pdbfixer \
    -y --solver=libmamba

./third_party/gmx_env/bin/pip install -r requirements.txt
```

驗證：
```bash
./third_party/gmx_env/bin/gmx --version | head -1     # GROMACS 2025.x
./third_party/gmx_env/bin/MMPBSA.py --version          # AmberTools MMPBSA
```

### 2. Tool A envs (Boltz / ColabFold / Chai-1)

```bash
# Boltz-2
conda create -p ./third_party/boltz_env python=3.10 -y --solver=libmamba
./third_party/boltz_env/bin/pip install boltz cuequivariance_torch cuequivariance-ops-torch-cu12
# 重要：Boltz 預設拉的 torch 對 CUDA 13；現場 driver 為 12.8，需手動換 wheel：
./third_party/boltz_env/bin/pip install --force-reinstall \
    torch --index-url https://download.pytorch.org/whl/cu128

# ColabFold (server mode — 不需本地 MMseqs2 資料庫)
conda create -p ./third_party/colabfold_env python=3.10 -y --solver=libmamba
./third_party/colabfold_env/bin/pip install 'colabfold[alphafold]'

# Chai-1
conda create -p ./third_party/chai_env python=3.10 -y --solver=libmamba
./third_party/chai_env/bin/pip install chai_lab
```

### 3. 模型權重路徑（可選）

預設 Boltz 把權重存到 `~/.boltz`（~5 GB）。若想留在 repo 內：
```bash
export BOLTZ_CACHE=$(pwd)/third_party/models/boltz
```

ColabFold 用 `--use-msa-server` 跑遠端 MSA，無需本地資料庫。Chai-1 同理（`--use-msa-server`）。

### 4. 受體準備（一次性）

```bash
# 下載 PDB 1J8E
mkdir -p dataset/lrp1_cr2/raw
curl -sL https://files.rcsb.org/download/1J8E.pdb -o dataset/lrp1_cr2/raw/1J8E.pdb

# 清洗 + 最小化 (Amber14 + OBC2 implicit solvent, ~10 秒於 A100)
./third_party/gmx_env/bin/python db/build_lrp1_receptor.py
```

成功會產出 `dataset/lrp1_cr2/receptor.pdb` + `metadata.json`，heavy-atom RMSD < 2 Å。

### 5. 預跑 demo 胜肽 cache（一次性）

三條 PAI-1 mimicking peptides（69-80, 76-88, 69-88）會被預跑進 `dataset/peptide_pai1/precomputed/`，之後 demo 時是 instant lookup。

```bash
export BOLTZ_CACHE=$(pwd)/third_party/models/boltz
./third_party/gmx_env/bin/python -c "
import sys; sys.path.insert(0, '.')
from engine.analyzer import ProteinOptimizer
opt = ProteinOptimizer()
for seq in ['KGMAPALRHLYK', 'RHLYKELMGPWNK', 'KGMAPALRHLYKELMGPWNK']:
    r = opt.run_pipeline(seq, tool_a='boltz2', tool_b='mmgbsa')
    opt.save_cache(seq, r)
    print(seq, r['binding_energy'])
"
```

每條約 3.5 分鐘（Boltz ~2 min + MM/GBSA ~1.5 min on A100）。

---

## Run

```bash
# 啟動後端
./third_party/gmx_env/bin/python api/server.py
```

API listens on `http://0.0.0.0:7860`. Frontend 走 Vite dev server（`cd frontend && npm run dev`）或 build 後由 Flask 服務。

### Pipeline behaviour
- 輸入是三條 demo 序列之一 → 直接從 `precomputed/` 回傳，毫秒內完成
- 其他序列 → 即時跑 Tool A → Tool B（A100 上約 3–8 分鐘，依序列長度）
- Tool A 預設 Boltz-2，可由 `tool_a` 參數切到 `colabfold` 或 `chai1`

### API

| Method | Path | 說明 |
|--------|------|------|
| POST | /api/jobs | 建立分析任務（body: `{sequence, tool_a?, tool_b?}`），回 `{job_id}` |
| GET  | /api/jobs/\<job_id\> | Polling 任務狀態與結果 |
| GET  | /api/uniprot/demo | UniProt demo targets |

---

## Notes

- **Tool B 為 single-point MM/GBSA**：sander GB minimization (igb=5, 2000 steps) 後 MMPBSA.py 取一格平均。Year-1 baseline，Year-2 會改全 MD ensemble。
- **Stefansson 校準**：`engine/analyzer.calculate_rho` 拿候選變體的 ΔΔG 對 Stefansson 1998 的 K69A/K80A/K88A Kd_fold 做 Spearman ρ。
- **受體選擇**：CR7（PDB 1J8E）為最佳解析度的 LRP1 CR 單一 domain。Year-2 升級到 CR3-CR10 tandem 或 AF2 模型。
