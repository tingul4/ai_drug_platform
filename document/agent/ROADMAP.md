# AIM3 胜肽藥物優化平台 ROADMAP

子計畫三：候選藥物之生成式設計與物理導向優化。標靶為 PAI-1 / LRP1 軸（胰臟癌）。
本文件統整目前實作狀態與下一步規劃，搭配 `README.md`（環境建置）與 `engine/analyzer.py`（pipeline 核心）閱讀。

---

## 1. 平台定位

針對 **PAI-1 mimicking peptides 競爭結合 LRP1 CR cluster II** 的設計問題，建立可信的計算流水線：

1. **3D 結構預測 (Tool A)** — Boltz-2 / ColabFold / Chai-1 三選一，輸出 peptide–receptor 複合物 PDB
2. **能量分析 (Tool B)** — AmberTools MM/GBSA（igb=5），輸出 ΔG_GB（kcal/mol）
3. **可靠度校準 (ρ)** — 模擬 ΔΔG vs Stefansson 1998 實驗 Kd_fold 的 Spearman 相關
4. **多目標排序 (MPO)** — Pareto ranking 跨結合能、免疫原性、聚集傾向

設計目標是 Year-1 baseline 跑通端到端，Stefansson 三筆實驗熱點（K69A / K80A / K88A）作為內部對照。

---

## 2. 系統架構

### 2.1 流水線（已串接，端到端可跑）

```
peptide sequence
    │
    ▼
[Tool A]  Boltz-2 (預設) / ColabFold / Chai-1
   - 受體序列固定為 LRP1 CR7（PDB 1J8E, 44 aa）
   - 走 ColabFold MSA server，無需本地 DB
   - 輸出：複合物 PDB（chain A = peptide, chain B = receptor）
    │
    ▼
[Tool B]  AmberTools MM/GBSA
   - tleap → ff14SB topology
   - sander GB minimization (igb=5, 2000 steps) — 解 Tool A 預測 pose 的微 clash
   - cpptraj → 1-frame netcdf
   - MMPBSA.py with `&gb` namelist (igb=5, saltcon=0.150)
   - 輸出：ΔG_GB (kcal/mol)
    │
    ▼
[Aggregation]
   - calculate_rho()：對候選變體做 Stefansson 校準
   - load_candidates() / run_alanine_scan()：目前為 stub（見 §3.2）
   - 結果寫入 cache 或回傳 JSON
```

### 2.2 環境（Multi-conda-env）

| Env | 內容 | 為什麼獨立 |
|---|---|---|
| `third_party/gmx_env` | GROMACS, AmberTools, gmx_MMPBSA, OpenMM, PDBFixer, Flask, biopython | 主 env，跑 Flask + Tool B |
| `third_party/boltz_env` | Boltz-2 + cuequivariance + torch (cu128) | numpy / biopython pin 與其他 Tool A 衝突 |
| `third_party/colabfold_env` | colabfold[alphafold] | 同上 |
| `third_party/chai_env` | chai_lab | 同上 |

`engine/analyzer.py` 透過 subprocess 跨 env 呼叫各 CLI。詳見 `README.md` 的 Environment setup。

### 2.3 後端 API（`api/server.py`）

| Method | Path | 功能 |
|---|---|---|
| POST | `/api/jobs` | body `{sequence, tool_a?, tool_b?}` → `{job_id}`，202 |
| GET | `/api/jobs/<id>` | polling status / progress / results |
| GET | `/api/uniprot/demo` | UniProt 示範 targets |

非同步 worker 執行 `optimizer.run_pipeline()`，cache 命中即時、未命中走 Tool A → Tool B 即時跑。

### 2.4 前端（`frontend/`）

React + Vite + Tailwind。`AnalysisPage.tsx` 提交序列、`ResultsPage.tsx` 輪詢結果並渲染 ρ、Pareto 表、Alanine scan。

---

## 3. 開發狀態

### 3.1 已完成 (Done)

#### Pipeline 核心
- [x] **Tool A 三條路全部接上**：Boltz-2 / ColabFold / Chai-1，均可由 `tool_a` 參數切換
- [x] **Tool B 實接 MM/GBSA**：tleap → sander minimization → MMPBSA.py，輸出物理合理的 ΔG_GB（demo peptides 落在 -46 ~ -69 kcal/mol）
- [x] **LRP1 CR7 受體準備**：`db/build_lrp1_receptor.py` 跑通（PDBFixer + OpenMM Amber14 + OBC2，heavy-atom RMSD = 0.683 Å）
- [x] **Stefansson ρ 校準邏輯**：`calculate_rho()` 接 K69A/K80A/K88A 對照
- [x] **Demo cache 機制**：`dataset/peptide_pai1/precomputed/<seq>.json`，命中直接回傳；三條 PAI-1 demo peptides 已預跑入庫
  - 69-80 (`KGMAPALRHLYK`)：ΔG_GB = -47.64 kcal/mol
  - 76-88 (`RHLYKELMGPWNK`)：ΔG_GB = -69.08 kcal/mol（最強，符合 Stefansson 主結合區）
  - 69-88 (`KGMAPALRHLYKELMGPWNK`)：ΔG_GB = -46.29 kcal/mol
- [x] **Live 路徑驗證**：任意輸入端到端跑通，11-mer 測試 243s 完成

#### 環境與基礎建設
- [x] **Multi-conda-env 架構**：四個獨立 env 解決 numpy/biopython pin 衝突
- [x] **GROMACS 2025.4 + AmberTools + gmx_MMPBSA**：bioconda 整包安裝
- [x] **Boltz torch CUDA 修正**：driver 12.8 → cu128 wheel（解 cu130 預設 wheel 的 driver too old 錯誤）
- [x] **README 環境建置文件**：完整重建步驟（含模型權重路徑、ColabFold MSA server 模式）
- [x] **Legacy SKEMPI 資產移除**：`benchmark/`、`candidates/`、`modifiable_sites/` 全清掉
- [x] **`runs/`、`third_party/`、`.venv/` 加入 .gitignore**

#### 後端 / 前端
- [x] **Flask Job Queue**：thread-pool worker、in-memory job table、progress 回報
- [x] **React + Vite 前端**：取代舊版 vanilla HTML/JS（`frontend/src/`）
- [x] **API 端到端串通**：`/api/jobs` 接 React 前端，cache hit 0.00s、live 走 polling

#### 周邊
- [x] **MHCflurry 2.0 + 台灣漢人 HLA-I panel**：免疫原性 scoring 模型整合（`engine/mhcflurry_immuno.py`）
- [x] **UniProt REST client + cache**（`engine/uniprot.py`）

### 3.2 已知限制（短期 backlog）

下列功能介面已預留，但目前回傳 stub 或硬編資料；需要在 `engine/analyzer.py` 替換實作：

- [ ] **`load_candidates()` 為 hardcoded 三筆**：固定回 K69A/K80A/K88A，未依輸入 peptide 動態產生
- [ ] **`run_alanine_scan()` 為 stub**：每位置回固定 ΔΔG=1.0、kd_fold=2.0，未實際 sweep
- [ ] **MPO 第二、三維度**：免疫風險、聚集傾向欄位在 candidates 為硬編字串
- [ ] **錯誤處理粗糙**：Boltz / sander / MMPBSA 任一段失敗，job status 變 error 但訊息不夠對使用者友善；無 GPU 排程或 timeout

---

## 4. 短期路線（Year-1 收尾，估計 1–2 週）

按優先順序：

### 4.1 前端端到端確認
- 啟動 Vite dev、實機跑三條 demo + 一條任意輸入
- 對齊 API JSON schema 與 `ResultsPage.tsx` 預期欄位
- 補 fallback：Boltz 失敗時前端顯示錯誤訊息（不只是空白頁）

### 4.2 `load_candidates` 改為 peptide-aware
- 依輸入序列實際殘基產生候選清單（K → A 替代僅在序列含 K 的位置）
- 保留快速 stub 路徑（不立即跑 Tool A/B）；ΔΔG 可先用 Statistical model 估
- 移除 hardcoded Stefansson rows，改由 `calculate_rho` 自己對應全域編號

### 4.3 Robustness
- subprocess 加 timeout 與清楚錯誤分類（Tool A fail / Tool B fail / parse fail）
- GPU 排程：避免兩個 Boltz job 同時搶 GPU 0；`CUDA_VISIBLE_DEVICES` 輪詢
- Job 超時自動標 `error` 並清理 `runs/<id>/`
- `/api/jobs/<id>` 增 `error_kind` 欄位給前端

### 4.4 真 Alanine scan（背景批次）
- 對每個位置跑單突變 Tool A + Tool B 取得真 ΔΔG
- 12-aa peptide 約 42 min；20-aa 約 70 min
- 結果存進 `precomputed/<seq>__ala.json` 並由 `run_alanine_scan` 讀取
- 三條 demo 全跑約 3 hr，可一次背景跑完

---

## 5. 中期路線（Year-2，~半年）

### 5.1 Pipeline 升級
- **真 candidate sweep**：每位置 19 種突變的 Tool A+B（成本：N × 19 × ~3.5 min；20-aa 一次 ~22 hr）；結合 Pareto ranking 做真 MPO
- **MD ensemble 取代 single-point**：sander short MD（~1 ns）→ 多 frame MMPBSA，提升 ΔG 統計穩定度
- **PSOS attention 加速**：長序列推理（如 70+ aa peptide variants）效能優化
- **多 Tool A 共識預測**：同時跑 Boltz / ColabFold / Chai 三者，取交集 / RMSD 過濾不可信 pose

### 5.2 受體升級
- 從 CR7 單一 domain 升級到 **CR3-CR10 tandem**（PAI-1 真實接觸的 LRP1 cluster II 完整範圍）
- 來源選項：AF2-multimer 預測、或拼接多個 CR PDB
- 重跑三條 demo cache 對齊新受體

### 5.3 MPO 完整化
- **免疫原性**：MHCflurry 2.0 串進 candidate generation 流程，每變體實際打分
- **聚集傾向**：AGGRESCAN / TANGO 串進序列分析
- **Pareto front 視覺化**：前端加 3D scatter（ΔG × immuno × agg）

---

## 6. 長期路線（Year-3+）

- **Year-3**：整合子計畫四 PBPK/QSP 臨床回饋迴路
- **Year-4**：端到端 SaaS 平台部署與多標靶適配（PAI-1/LRP1 之外擴到其他癌症 PPI）
- **遠期**：active learning loop — 計算結果回灌實驗驗證、再校準 ΔG 預測

---

## 7. 文件版本

- 2026-04-25：本次大改版。Section 3 全面更新（Tool A/B 真接、demo cache、conda 環境）。短期路線改為具體 4.1–4.4，原本 Year-2/3/4 框架保留為 §5–6。
- 2026-04-24：合併取代 `Aim3_POC_implementation.md`、`AIM3_成果條列_v2.md`、`DEMO_AIM3.md`、`peptide_pipeline_plan.md`、`roadmap.md`、`v1_analysis_summary.md`。
