# 子計畫三 — 實作 Roadmap

依據 `document/human/工程CM03-結合-unify-20260201-v35.pdf`（計畫書）第 25-28 頁
「子計畫 3：AI 藥物設計與最佳化 — 以 PAI-1 之胜肽先導分子為主案例」，
與 `document/agent/peptide_pipeline_plan.md` 定義的 Year-1 Baseline Pipeline。

對應標的：**三條 PAI-1 mimicking peptides（69-80, 76-88, 69-88）→ LRP1 CR cluster II (PDB 1J8E)**。

---

## 總覽

| 階段 | 負責範圍 | 完成度 |
|---|---|---|
| Year-1：Baseline Pipeline & 資料基礎 | S1-S8（見下） | 🟡 進行中 |
| Year-2：AI-Physics Hybrid + PSOS 加速 | 瓶頸處導入 PSOS / Deep mutation scan | ❌ 未開始 |
| Year-3：多參數整合 + 臨床回饋 | 子計畫 4 PBPK/QSP 接口 | ❌ 未開始 |
| Year-4：端到端部署 + 平台化 | 整合 Demo + SaaS | ❌ 未開始 |

---

## Year-1 Baseline — 8 步驟 MVP

對齊計畫書圖四「Four-Year AI Drug Design & Optimization Plan」第一年橫條
（Data Harmonization → Structural Modeling → Structural Modeling Baseline → Basic Scoring & Screening → Initial Candidates）。

| # | 步驟 | 產出 | 驗證條件 | 狀態 |
|---|---|---|---|---|
| S1 | Data Harmonization | `dataset/peptide_pai1/` FASTA + JSON metadata（3 條胜肽）+ `db/build_peptide_dataset.py` | schema 驗證通過、保留 FITC-Ahx/MW/purity/batch metadata | ✅ MW Δ=0.00 通過 |
| S2 | Target 結構準備 | `dataset/lrp1_cr2/receptor.pdb`（PDB 1J8E + PDBFixer + Amber14/OBC2 min） | 最小化收斂；RMSD to crystal < 2 Å | ✅ RMSD 0.551 Å |
| S3 | Complex Modeling | AF-Multimer：3 peptide × LRP1 CR2；pLDDT + 介面殘基匯出 | interface pLDDT > 70；界面殘基 ∩ 文獻 PAI-1–LRP1 接觸 > 50% | ❌ |
| S4 | Peptide Docking 補強 | HPEPDOCK / ADCP 獨立 docking；PLIF 指紋；Alanine scan | 與 `SKEMPIModel` 預測 Spearman > 0.6 | ❌ |
| S5 | Modifiable Sites List | 固定 schema JSON：`residue_index, allowed_mutations, interface_importance_score, interaction_type, evidence_source` | schema conform；`residue_index ∈ [1, N]` | ❌ |
| S6 | Immunogenicity | MHCflurry 2.0（台灣 HLA-I panel）主 + IEDB-MHC-II PSSM 副 | 三條 peptide 的 MHC-I %rank 分布可畫；`%rank<2` 計數標註 | ✅ combined = 0.7·I + 0.3·II |
| S7 | ADMET + Tox + Synth | ADMETlab 2.0 REST + Tox21 分類 + AiZynthFinder 合成路徑 | 三條 peptide 拿到完整 ADMET panel + 合成評分 | ❌ |
| S8 | J(x) Baseline | `engine/analyzer.py` 擴充；統一權重 `w_i = 0.25`；Penalty 硬條件：any window MHC-I %rank<2 or Tox21 flag | J(x) 可排序、有解釋；回歸測試通過 | ❌ |

---

## 現況沿用（無需重寫）

以下 v1 元件仍是 Year-1 baseline 的核心，**不需要重做**：

- `SKEMPIModel` / `predict_ddG` / `predict_koff`（`engine/analyzer.py:204-`）
- `alanine_scan` — S4 依賴
- `generate_variants` / `score_variants` — Year-2 Deep mutation scan 的起點
- NSGA-II Pareto（`_fast_nds`） — J(x) 多目標排序
- SQLite session 持久化（`db/skempi.db`）
- Flask REST API + 前端（`api/server.py`, `frontend/index.html`）
- UniProt 快取（`engine/uniprot.py`） — S2 LRP1 抓取會沿用
- IEDB-MHC-II PSSM（`engine/analyzer.py:134-200`） — **降為副指標**，不移除

---

## 現況調整 / 重做

- `calc_immunogenicity`（heuristic）**移除**，已被 PSSM 取代；且計畫書明確用 NetMHCpan/MHCflurry。
- `score_immunogenicity` 回傳值：scalar → dict `{mhc_i_risk, mhc_ii_risk, combined, top_hits}`（詳 `peptide_pipeline_plan.md §4.5`）。
- Pareto 三目標（ΔΔG / immuno / aggregation）→ 擴為 J(x) 四項
  `J = w₁·f_bind + w₂·f_admet + w₃·f_synth − w₄·Penalty`（計畫書 p.28）。
- `uniprot.py` 抓取範圍：除 PAI-1 (P05121) 外，新增 LRP1 (P98157) 的 CR cluster II 切片。

---

## ❌ 明確 Out-of-Scope（不做）

下列項目在舊版 roadmap（`document/Aim3_doc_20260119.pdf` 對應）中曾列為目標，
但依計畫書 CM03 v35 第 25-28 頁的工作分工，**不屬於子計畫 3 的 Year-1 工作**，
或已被更精準的工具取代：

### 小分子分支

- KRAS G12D 小分子 POC（`engine/smallmol.py`、`dataset/kras_g12d_poc/`）
  → 歷史 POC，保留程式碼但**不納入 Year-1 MVP**（見 `Aim3_POC_implementation.md` 頂部狀態註）。
- Uni-Mol / Uni-Mol+ 小分子 3D 表徵。
- RDKit ETKDGv3 + Meeko + AutoDock Vina（小分子）docking 流程。

### 屬於子計畫 2（Council of Agents）的工作

- DrugBank / PubMed / ChEMBL 自動化檢索 agent（A3）。
- Open Targets GraphQL API / Bayesian 信心評分（A2）。
- Ensembl / dbSNP / GTEx / DisGeNET / OMIM / STRING / WGCNA 工具鏈（A1）。
- Main Orchestrator / ReAct / ToolOrchestra / Chain-of-Evidence。
- LLM Markdown 敘事報告、Supervised Reinforcement Learning (SRL)。

### 屬於子計畫 4（PBPK/QSP 虛擬臨床）的工作

- PK-Sim / MoBi / PBPK / QSP 模擬。
- Causal Forest 族群分層。
- Digital Twin / 虛擬臨床試驗。

### Year-2 以後才做（Year-1 只做 profiling）

- **PSOS (Prefix-Scan Online Softmax) 演算法**：計畫書 p.28「本年度結束時…Profiling 結果將在第二年導入 PSOS 前形成明確的瓶頸定位與量化目標」。
- Deep Mutation Scanning（全序列突變掃描）。
- AI-Physics Hybrid Optimization。
- Finer Structure Scoring（候選–標靶複合體局部鬆弛）。
- Interpretability Module（Attention map、SHAP、3D 歸因圖譜）。

### 其他被取代 / 非必要項目

| 舊 roadmap 項目 | 處理方式 | 替代方案 / 原因 |
|---|---|---|
| ESM-2 殘基級嵌入 | ⛔ 不做 | 未列於計畫書子計畫 3 工具鏈；複合結構由 AF-Multimer 提供更直接的特徵 |
| FoldX | ⛔ 不做 | 授權卡商用；SKEMPIModel 已涵蓋 ΔΔG baseline，Year-2 升級走 AF-Multimer + Rosetta |
| τ_RAMD / AI-MD / OpenMM + TorchMD-NET | ⛔ Year-1 不做 | 計畫書 Year-1 不要求 MD；結合動力學以 SKEMPI koff_ratio 估算已足夠 baseline |
| AGGRESCAN 3D / CamSol 官方工具 / NetSurfP-3.0 / SCM | ⛔ 不做 | 現有 heuristic 已涵蓋；真正的 aggregation / solubility 評估由 ADMETlab 2.0 統一提供 |
| FcRn 動態結合模擬 | ⛔ 不做 | 抗體藥物 PK 範疇，屬子計畫 4 |
| AutoDock Vina（用於胜肽） | ⛔ 改 HPEPDOCK / ADCP | Vina 對 > 6-7 AA 的 flexible peptide 旋轉自由度爆炸、不準 |
| NetMHCpan-4.1（bundle 進 repo） | ⛔ 不 ship | 學術授權禁商用；若本地有授權則當 internal reference adapter 校準 MHCflurry threshold |
| MHC-I PSSM（自建） | ⛔ 不做 | MHC-I groove 封閉、anchor 位置耦合強，PSSM 在 MHC-I 準確度不夠；直接上 MHCflurry 2.0 |
| LLM Markdown 報告、PyMOL 視覺化 | ⛔ Year-1 不做 | 屬於 Year-4 Demo 階段 |
| 實驗閉環回饋（Aim 3C） | ⛔ Year-1 不做 | 計畫書寫在 Year-3「臨床回饋介面」 |

---

## 依賴與風險

- **LRP1 結構缺口**：沈延盛教授團隊目前無內部 PAI-1–LRP1 共結晶結構 → 採公開 PDB 1J8E（LRP1 CR cluster II）為 receptor。若後續團隊有 co-crystal，S3 AF-Multimer 結果需重新校準。
- **MHCflurry 2.0 模型權重**：一次性下載 ~500 MB；repo 用 uv 管理 `.venv/`，安裝命令見 `peptide_pipeline_plan.md §5 Q2`。
- **HLA panel**：採台灣漢人高頻 HLA-A/B/C allele（與計畫書 p.8「情境化校準」論述對齊），panel 列表見 `peptide_pipeline_plan.md §5 Q3`。
- **ADMETlab 2.0**：為第三方 REST API，需評估 rate limit 與離線 fallback；Tox21 可用 DeepChem 本地模型。
- **AiZynthFinder**：需下載 policy network 模型（~1 GB），若磁碟緊張可先跳過合成評分、`f_synth` 預設為 0.5。
