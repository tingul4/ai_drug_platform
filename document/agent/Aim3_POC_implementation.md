# Aim3 POC 實作紀錄：KRAS G12D 胰臟癌候選分子篩選 Pipeline

> **⚠️ 狀態：歷史紀錄（Historical / Out-of-scope for Year-1）**
>
> 此份文件記錄 2026-Q1 基於 `Aim3_doc_20260119.pdf` 所做的**小分子 KRAS G12D POC**。
> 在 2026-04 拿到實際的 PAI-1 mimicking peptides（`Peptide info.zip`）後，
> 子計畫三第一年的主線已根據 `工程CM03-結合-unify-20260201-v35.pdf`
> 與 `peptide_pipeline_plan.md` 重新定位為 **PAI-1 胜肽 → LRP1 受體**路線。
>
> **本 POC 之角色：** 保留為「小分子旁路」的可行性證明與 `engine/smallmol.py` 程式碼的設計依據，
> **不納入第一年 Year-1 MVP 交付範圍**。若未來有 KRAS 或其他小分子標靶需求再啟用。
>
> 現行第一年 roadmap 請見：`document/agent/roadmap.md` 與 `document/agent/peptide_pipeline_plan.md`。

---

> 對應文件：`document/Aim3_doc_20260119.pdf`  
> 目標：建立一版「功能完整、不複雜」的端到端 pipeline，重現 PDF 預期輸出的 Summary Table 與 Pareto 排序。

---

## 1. 總覽：PDF 預期輸入輸出 ↔ 本 POC 實作對照

| PDF 項目 | 本 POC 對應檔案 / 工具 |
|---|---|
| 候選藥物 SMILES | `dataset/kras_g12d_poc/ligands/candidates.tsv`（10 個分子）|
| 標靶蛋白 FASTA / PDB | `structure/KRAS_P01116.fasta` + `structure/7RPZ.pdb` |
| pIC50 / pKd 預測 | AutoDock Vina 1.2.3 → `-ΔG/(2.303RT)` proxy |
| koff 動力學代理 | Descriptor-based proxy（HBD、ar rings、rot bonds、MW）|
| LogP / hERG / CYP3A4 | ADMET-AI（Graph Attention Network）|
| Pareto 多目標排序 | PyMOO NonDominatedSorting + NSGA-II |
| Summary Table | `results/summary_table.tsv` |
| Pareto Front 輸出 | `results/pareto_front.tsv` + `pareto_plot.png` |

**候選分子（10 個 POC）**：Sotorasib、Adagrasib、MRTX1133、Divarasib、BI-2852、Osimertinib、Gefitinib、Erlotinib、Imatinib、Dasatinib。前 5 個為 KRAS 抑制劑 parent/變體，後 5 個為已上市 kinase 藥物作為對照。

---

## 2. 環境建立

**做了什麼**：用 `uv` 建立 `.venv`（Python 3.10），以 `requirements.txt` 管理依賴。

**為什麼**：
- `uv` 比 `pip/conda` 快一個數量級，且可重現。
- Python 3.10 是 chemprop（ADMET-AI 的底層）相容性最穩定的版本；3.13 有 ABI 問題。
- 依賴清單見 `requirements.txt`。關鍵 pin：`numpy>=1.26,<2.0`（vina 1.2.x 對 numpy 2.x ABI 不相容）。

**產出**：`.venv/`、`requirements.txt`

---

## 3. 標靶準備：7RPZ + KRAS FASTA

**做了什麼**：`curl` 抓 RCSB PDB `7RPZ`（KRAS G12D + MRTX1133 共結晶）與 UniProt `P01116` FASTA。

**為什麼選 7RPZ**：
- G12D 是 PDAC 最常見的 KRAS 突變（約 35% 胰臟癌）。
- 7RPZ 提供 G12D ATP-site 抑制劑 MRTX1133 的已知結合姿態 → 可直接當 pocket anchor，不用再做 pocket detection。
- 解析度 1.38 Å（高品質），適合 docking。

**產出**：`structure/7RPZ.pdb`、`structure/KRAS_P01116.fasta`

---

## 4. 候選 SMILES 抓取與 3D 結構生成

**做了什麼**：
1. `ligands/fetch_smiles.sh`：用 PubChem REST API（`/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/TXT`）抓 canonical SMILES，避免手打出錯。
2. `ligands/prepare_ligands.py`：對每個 SMILES 做 RDKit `ETKDGv3` embed 生成 3D + UFF 最佳化，寫 SDF 與 Meeko PDBQT。

**為什麼這樣做**：
- **PubChem CID 而非手動 SMILES**：確保 SMILES 正確，失誤風險低。
- **ETKDGv3**：比舊版 ETKDG 更穩；對含立體中心/環系分子成功率更高。
- **UFF 簡單 min**：POC 夠用；要更精準可換 MMFF94 或 xTB。
- **Meeko PDBQT**：Vina 需要 PDBQT。Meeko 自動加 Gasteiger 電荷、合併非極性 H、生成 torsion tree，一行搞定。

**產出**：`ligands/candidates.sdf`、`ligands/pdbqt/CAND_0{1..10}.pdbqt`（10/10 成功）

---

## 5. Docking：Vina 1.2.3 CLI + 原生 ligand anchor

**做了什麼**：`docking/run_docking.py`
1. 解析 7RPZ.pdb 中的 MRTX1133（resname `6IC`）座標，算 centroid 當 pocket center。
2. `mk_prepare_receptor.py` 把只保留 protein 的 PDB 轉 PDBQT。
3. `/usr/bin/vina`（CLI，v1.2.3）對每個 ligand 跑 docking，exhaustiveness=16、n_poses=10。
4. 取最佳 affinity（kcal/mol），用 `pKd = -ΔG / (2.303·RT)`（298 K 下 = ΔG / 1.3633）換算成 pKd proxy；pIC50 POC 假設等於 pKd（競爭性結合）。

**為什麼這樣做（重要 debugging 紀錄）**：
- 一開始用 `vina` Python binding（1.2.5 / 1.2.7）失敗，錯誤「Affinity map for atom type Cl_H / N_D is not present」。原因：新版 Vina 在 `dock()` 時會把 ligand 的 Cl/N 再分類為 `Cl_H`（halogen bond donor）、`N_D`（H-bond donor），但 `compute_vina_maps()` 只根據 PDBQT 原始型別建 map，造成 map 不對應。
- **解法**：改用系統的 `/usr/bin/vina` CLI（v1.2.3，halogen bonding 前的版本），atom type 處理較單純；Meeko 預設的 AD4 type 直接對得上。這也記錄在腳本開頭 docstring 中，避免未來踩同樣坑。
- **Pocket center 用共結晶 ligand centroid**：不需 pocket detection，最準且最快。

**產出**：`docking/docking_scores.tsv`、`docking/CAND_0{1..10}_poses.pdbqt`、`docking/pocket_center.json`、`docking/receptor.pdbqt`

**結果摘要**（最佳 5 名）：
| ID | 分子 | ΔG (kcal/mol) | pKd proxy |
|---|---|---:|---:|
| CAND_09 | Imatinib | −9.87 | 7.24 |
| CAND_10 | Dasatinib | −9.58 | 7.02 |
| CAND_06 | Osimertinib | −9.16 | 6.72 |
| CAND_02 | Adagrasib | −8.84 | 6.48 |
| CAND_08 | Erlotinib | −8.43 | 6.18 |

> Imatinib/Dasatinib 排第一第二並非 KRAS 原生活性高，而是它們是多靶點 kinase 抑制劑，KRAS ATP site 對類似骨架有一定親和力；這正好測試 pipeline 能「排序」的能力。

---

## 6. ADMET-AI：一次預測 LogP / hERG / CYP3A4 等 80+ 指標

**做了什麼**：`admet/run_admet.py` 呼叫 `admet_ai.ADMETModel()` 對 10 個 SMILES 做 batch 預測。

**為什麼選 ADMET-AI**：
- 核心技術：多任務 **Graph Attention Network (GAT)**，同時學 80+ 指標間的關聯，不只看單一性質。
- 一個 API 涵蓋 Summary Table 的三個欄位（LogP、hERG、CYP3A4）以及 BBB、AMES、DILI、ClinTox 等 bonus。
- PDF 明確列為「強烈推薦」。

**Debugging 紀錄**：PyTorch ≥2.6 把 `torch.load` 預設 `weights_only=True`，但 ADMET-AI/chemprop 的 checkpoint 裡有 `argparse.Namespace`、`numpy.core.multiarray._reconstruct` 等非白名單 globals，會 `UnpicklingError`。
**解法**：在 import `admet_ai` 前 monkey-patch `torch.load` 強制 `weights_only=False`（ADMET-AI 官方 checkpoint 可信）。

**產出**：`admet/admet_scores.tsv`（100 欄全指標）、`admet/admet_summary.tsv`（關鍵欄位）

---

## 7. koff proxy：描述子啟發式代理，不跑 MD

**做了什麼**：`admet/koff_proxy.py`，公式：
```
proxy = sigmoid(0.40·HBD + 0.25·aromatic_rings − 0.15·rot_bonds − 0.003·(MW − 450))
≥ 0.60 → Low (Stable)   (slow koff, long residence time — preferred)
0.40 ~ 0.59 → Medium
< 0.40 → High (Fast)
```

**為什麼不跑真 τRAMD**：
- 真 τRAMD 每個 ligand 需要 10–20 個 replicate × 數分鐘 GPU，10 個分子需數小時，不符合「POC 先不複雜」的要求。
- Residence time 文獻（Schmidtke 2011, Copeland 2016）指出 slow off-rate 與 HBD、aromatic stacking、剛性骨架（低 rot_bonds）正相關 → 以這些描述子做線性啟發式代理是可接受的 POC 近似。
- **明確標示「不是真 kinetic 量測」**；未來升級路徑：改呼叫 `OpenMM + τRAMD` pipeline。

**產出**：`admet/koff_proxy.tsv`

---

## 8. Summary Table 組合

**做了什麼**：`results/build_summary.py` 把三個 TSV join 起來，欄位對齊 PDF：`ID | pIC50 | pKd | koff | LogP | hERG 風險 | CYP3A4 抑制`。

**風險標籤規則**：
- hERG：機率 ≥ 0.5 → `High`，否則 `Low`
- CYP3A4：機率 ≥ 0.5 → `Yes`，否則 `No`

**產出**：`results/summary_table.tsv`（10 列，對應 PDF 範例表格格式）

---

## 9. PyMOO NSGA-II 多目標排序

**做了什麼**：`results/run_nsga2.py` 分兩階段：

### Stage A：真實候選的 Pareto 排序（主要成果）

用 `pymoo.util.nds.non_dominated_sorting.NonDominatedSorting` 直接對 10 個候選做非支配排序，**6 目標全部 minimize**：

| 目標 | 意義 |
|---|---|
| `-pKd` | 親和力越高越好（取負） |
| `-pIC50` | 同上 |
| `-koff_score` | koff proxy 分數越高 koff 越慢（取負） |
| `\|LogP − 2.5\|` | drug-like 甜蜜點 |
| `hERG` | 心毒性機率，越低越好 |
| `CYP3A4` | 代謝抑制機率，越低越好 |

### Stage B：NSGA-II 在理想性質空間的示範搜索

對 `(LogP, koff_score)` 2D 連續空間用 Ridge 代理 6 目標，NSGA-II 80 代 × pop=80，繪出「理想的 Pareto Front 形狀」供未來分子生成參考。

**為什麼選 NSGA-II**：
- 6 目標在 NSGA-II 擅長區間（2–10 目標）
- PyMOO 內建、成熟、參考實作豐富
- 若日後擴充到 >10 目標再換 NSGA-III / MOEA/D

**Stage A 結果（Pareto Front 5/10）**：

| ID | 分子 | pKd | LogP | hERG | CYP3A4 | koff_score | 被選原因 |
|---|---|---:|---:|---:|---:|---:|---|
| CAND_09 | Imatinib | 7.24 | 4.59 | 0.97 | 0.86 | 0.65 | 親和力最高 |
| CAND_10 | Dasatinib | 7.02 | 3.31 | 0.93 | 0.51 | 0.69 | 親和力 + koff 平衡 |
| CAND_08 | Erlotinib | 6.18 | 3.41 | 0.92 | 0.87 | 0.46 | LogP 近甜蜜點 |
| CAND_05 | BI-2852 | 5.73 | 3.11 | 0.43 | 0.15 | 0.61 | **hERG/CYP3A4 極低** |
| CAND_03 | MRTX1133 | 5.68 | 2.54 | 0.34 | 0.24 | 0.67 | **毒性最安全** |

**解讀**：Pareto Front 把「高親和力但有毒」（Imatinib）與「低毒但親和力中等」（MRTX1133）都保留，正是多目標最佳化的價值——**沒有單一最好，只有一組權衡最佳**，讓下游決策者依優先權挑選。

**產出**：`results/pareto_front.tsv`、`results/candidates_ranked.tsv`、`results/pareto_plot.png`、`results/nsga2_idealized_front.tsv`

---

## 10. 可視化：Pareto plot

`results/pareto_plot.png` 雙欄散佈圖：
- 左：pKd vs hERG（親和力 vs 心毒性）
- 右：pKd vs koff_score（親和力 vs 停留時間）
- 顏色 = Pareto rank（0 = 非支配）

快速看出「Pareto 前沿」上哪些分子值得進一步實驗驗證。

---

## 11. 對照 PDF Summary Table 欄位

| PDF 欄位 | 本 POC 欄位 | 檔案 |
|---|---|---|
| 候選編號 (ID) | `id` | summary_table.tsv |
| pIC50 (活性預測) | `pIC50` | summary_table.tsv |
| pKd (親和力) | `pKd` | summary_table.tsv |
| koff 傾向 (動力學) | `koff` (Low/Medium/High) | summary_table.tsv |
| LogP (理化) | `LogP` | summary_table.tsv |
| hERG 風險 (心毒性) | `hERG_risk` (Low/High) | summary_table.tsv |
| CYP3A4 抑制 (代謝風險) | `CYP3A4_inhibition` (Yes/No) | summary_table.tsv |

PDF 要求的 **全方位去風險評估矩陣**：`admet/admet_scores.tsv` 還有 BBB、AMES、DILI、ClinTox 等 bonus 欄位可擴充使用。

---

## 12. 端到端執行

```bash
source .venv/bin/activate
cd dataset/kras_g12d_poc

# 1. 下載結構
cd structure && curl -sL https://files.rcsb.org/download/7RPZ.pdb -o 7RPZ.pdb && \
  curl -sL https://rest.uniprot.org/uniprotkb/P01116.fasta -o KRAS_P01116.fasta

# 2. 候選分子
cd ../ligands && bash fetch_smiles.sh && python prepare_ligands.py

# 3. docking
cd ../docking && python run_docking.py

# 4. ADMET + koff proxy
cd ../admet && python run_admet.py && python koff_proxy.py

# 5. Summary + NSGA-II
cd ../results && python build_summary.py && python run_nsga2.py
```

## 13. 下一階段 TODO

對齊 PDF 的三個階段：

### 第一階段強化（Multi-modal representation）
- [ ] 以 **Uni-Mol** 生成 3D embedding，替代 RDKit ETKDG
- [ ] **ESM-2** 抽 KRAS 殘基級特徵，用於口袋/熱點識別
- [ ] **AlphaFold2 / ESMFold** 為缺結構標靶建模

### 第二階段強化（Kinetic & Selectivity）
- [ ] **真 τRAMD**（OpenMM + `kofflist`）取代 descriptor proxy
- [ ] 加入 **3D-GNN 注意力**（PyG）辨識關鍵 residue interactions
- [ ] **SHAP** 歸因分析：哪些子結構驅動了 hERG risk

### 第三階段強化（De-risking）
- [ ] 蛋白分支：整合 ESM-2 + **FoldX** 預測 Tm、**AGGRESCAN 3D** 聚集、**NetMHCpan** 免疫原性
- [ ] 數位風險標籤自動分級（Rule-based + SHAP）
- [ ] LLM 敘事報告（LLaMA-3 本地部署），產出 PDF 範例中的「AI 驅動先導化合物優化報告」

### 規模擴充
- [ ] 從 **ZINC20** / **GDB-17** 拉數萬個候選跑 virtual screening
- [ ] 用 **REINVENT 4.0** 做生成式分子設計，補足 Stage B 的理想分子
- [ ] 改用 **NSGA-III / MOEA/D** 若目標數 >10

### 資料庫擴充（若需訓練自家模型）
- [ ] 下載 **ChEMBL** / **BindingDB** 訓練客製親和力模型
- [ ] 下載 **MoleculeNet** 做 ADMET 模型 benchmarking

---

## 附錄 A：檔案樹

```
dataset/kras_g12d_poc/
├── structure/
│   ├── 7RPZ.pdb
│   └── KRAS_P01116.fasta
├── ligands/
│   ├── fetch_smiles.sh
│   ├── prepare_ligands.py
│   ├── candidates_raw.tsv
│   ├── candidates.tsv
│   ├── candidates.sdf
│   └── pdbqt/CAND_0{1..10}.pdbqt
├── docking/
│   ├── run_docking.py
│   ├── receptor.pdbqt
│   ├── pocket_center.json
│   ├── CAND_0{1..10}_poses.pdbqt
│   └── docking_scores.tsv
├── admet/
│   ├── run_admet.py
│   ├── koff_proxy.py
│   ├── admet_scores.tsv      (100 cols)
│   ├── admet_summary.tsv
│   └── koff_proxy.tsv
└── results/
    ├── build_summary.py
    ├── run_nsga2.py
    ├── summary_table.tsv     ← 對應 PDF Summary Table
    ├── candidates_ranked.tsv
    ├── pareto_front.tsv      ← 5/10 非支配分子
    ├── nsga2_idealized_front.tsv
    └── pareto_plot.png
```

## 附錄 B：關鍵版本

- Python 3.10.0 / uv 0.10.2
- RDKit 2026.03.1
- PyMOO 0.6.1.6
- AutoDock Vina 1.2.3 (CLI) + 1.2.5 (Python binding, unused)
- Meeko（mk_prepare_receptor.py / mk_prepare_ligand.py）
- ADMET-AI（chemprop backend）
- PyTorch 2.11.0
- NumPy 1.26.4（pin for vina ABI compatibility）
