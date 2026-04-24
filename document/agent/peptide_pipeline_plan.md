# Peptide Pipeline 第一年 Baseline 規劃

作為 `document/human/response.md` 中「第一年實作建議」的正式分析與確認稿。
對照來源：`document/human/工程CM03-結合-unify-20260201-v35.pdf`（以下稱計畫書）、
`Peptide info.zip`（GenScript 合成 COA）、`engine/analyzer.py` 現況。

---

## 1. 手上的三條胜肽（來自 `Peptide info.zip`）

| Product | Order ID | Sequence | Length | N-term Mod | MW | HPLC purity |
|---|---|---|---|---|---|---|
| 69-88 | U230SGAJG0_7  | `KGMAPALRHLYKELMGPWNK` | 20 AA | FITC-Ahx | 2843.36 | 98.2% |
| 69-80 | U230SGAJG0_9  | `KGMAPALRHLYK`         | 12 AA | FITC-Ahx | 1887.24 | 98.0% |
| 76-88 | U230SGAJG0_11 | `RHLYKELMGPWNK`        | 13 AA | FITC-Ahx | 2174.51 | 99.3% |

對應計畫書第 1 頁「PAI-1 mimicking peptides（6980、7688、6988）」。

---

## 2. ⚠️ 關鍵修正：對接標靶是 LRP1，不是 PAI-1

計畫書第 1 頁原文：

> 另也設計 PAI-1 mimicking peptides（共有三個候選胜肽 6980、7688、6988），
> 初步顯示可與 **LRP1 結合並被 PAI-1 競爭**。

這三條胜肽是**模仿 PAI-1 的 LRP1 結合片段**，目的是**競爭性佔住 LRP1 上的 PAI-1 結合位**。
所以：

- Docking / complex modeling 的 **receptor = LRP1**（不是 PAI-1）。
- `engine/uniprot.py` 原本抓 PAI-1 (SERPINE1, P05121) 是抓「被模仿者」，在 docking 階段要另外準備 LRP1 (UniProt P98157) 的 CR cluster 結構。
- N-term FITC-Ahx 是螢光探針（團隊做 LRP1 結合 assay 用），**電腦模擬請用裸序列**，但 metadata 保留修飾資訊與 batch ID。

---

## 3. 針對 `document/human/response.md` 5 步驟的逐項評估

| 步驟 | 原建議 | 計畫書 p.26-28 | 修正 / 建議 |
|---|---|---|---|
| S1 資料標準化 | FASTA | ✅「胜肽以 FASTA 為主，必要時搭配 PDB」 | 加入固定 schema：`modifiable_sites_list` 必含 `residue_index`, `allowed_mutations`, `interface_importance_score`, `interaction_type`, `evidence_source` |
| S2 結構建模 | AlphaFold / ESMFold + CHARMM | ✅「以 AlphaFold [23] 與 ESMFold [24] …CHARMM/CGenFF 最小化」 | **對 12-20 AA 裸胜肽單獨跑 AlphaFold 通常低 pLDDT 無意義**。改為 peptide + LRP1 CR domain 的 **complex modeling**（AlphaFold-Multimer 或 HADDOCK with interface restraints） |
| S3 Docking + PLIF + Ala scan | AutoDock Vina 半柔性 | 只寫「Docking / 能量評分」+ PLIF + Alanine scanning | **Vina 對 > 6-7 AA 的 flexible peptide 旋轉自由度太大、不準**。改 **HPEPDOCK 2.0** / **AutoDock CrankPep (ADCP)** / AF-Multimer 複合結構 + **Rosetta FlexPepDock** refine |
| S4 免疫原性 + ADMET | NetMHCpan-4.1 + ADMETlab 2.0 | ✅ 原文指定 NetMHCpan [30] + ADMETlab 2.0 [29] + Tox21 [32] + AiZynthFinder [28] | **計畫書還要求 Tox21（毒性分類補強）與 AiZynthFinder（合成可行性）**，repo 目前都沒有 |
| S5 目標函數 J(x) | w₁·f_bind + w₂·f_admet + w₃·f_synth − w₄·Penalty | ✅ 原文一致 | Penalty 硬條件至少含「明顯毒性警訊」與「免疫原性風險過高」。PSOS 第二年才導入，第一年不要動 |

---

## 4. MHC-I vs MHC-II 與工具選型

### 4.1 結論

> 對這批 12-20 AA 的短線性胜肽，應**以 MHC-I 為主、MHC-II 為副**。
> **主要工具選 MHCflurry 2.0**（Apache-2.0、SOTA 準確度、Python-native、可商用）。
> **不自建 MHC-I PSSM**。**不優先上 NetMHCpan**（雖計畫書點名，但授權會卡後續 SaaS/商用）。
> 現有 `IEDBImmunoModel`（MHC-II PSSM, `engine/analyzer.py:134-200`）降級為副指標保留。

### 4.2 與前一輪「生物藥去免疫化偏 MHC-II」的差異

前一次討論的是 ~150-450 AA 的抗體/蛋白藥，ADA 誘導經 CD4+ helper pathway，所以 MHC-II 重要。
這批胜肽不同：

- 長度 12-20 AA，APC 吸收後直接進 proteasome → **MHC-I (CD8+) 呈遞路徑為主**。
- 為 PAI-1 mimetic，序列本身來自人類自體蛋白，central tolerance 會抑制原生 epitope，
  但若後續做突變優化引入非 self 側鏈，**MHC-I 突破 tolerance 的風險反而是首要關注**。
- MHC-II 還需要看（CD4 幫助 ADA），但退到第二順位。

### 4.3 為什麼 MHC-I 不做 PSSM

- MHC-I binding groove 兩端封閉，**P2 與 P9 anchor 與中間殘基有強交互作用**，PSSM 的位置獨立假設在 MHC-I 會吃大的準確度虧。
- MHC-II groove 兩端開放，PSSM 堪用（這也是現有 IEDB-MHC-II PSSM 有基本效度的原因）。
- 這也是近十年 MHC-I 工具全面改 NN（NetMHCpan、MHCflurry）的主因。

### 4.4 工具比較

| 工具 | 型別 | License | MHC-I AUC | Python 整合 | 適用度 |
|---|---|---|---|---|---|
| NetMHCpan-4.1 | NN | 學術授權（禁商用） | ~0.93 | binary wrapper | ⚠️ 可當 reference，不當主力 |
| **MHCflurry 2.0** | NN (TF/Keras) | **Apache-2.0** | ~0.92 | `pip install mhcflurry` | ✅ **主要工具** |
| MHC-I PSSM（自建） | PSSM | 無 | ~0.75-0.82 | 現有程式可延伸 | ❌ 準確度不夠 |
| MixMHCpred | PSSM+ | 開源但禁商用 | ~0.90 | CLI | ⚠️ 商用卡關 |

### 4.5 落地到 `engine/analyzer.py`

現在 `score_immunogenicity` (`engine/analyzer.py:374`) 只回 scalar。改為：

```python
immuno = {
    "mhc_i_risk":   mhcflurry_predictor.score(seq),   # 主，per-allele % rank → risk
    "mhc_ii_risk":  self.iedb_pssm.score(seq),         # 副，保留現況
    "combined":     0.7 * mhc_i_risk + 0.3 * mhc_ii_risk,  # 給 NSGA-II
    "mhc_i_top_hits":  [...],
    "mhc_ii_top_hits": [...],  # 現有 PSSM 輸出
}
```

硬條件 Penalty：任何 9-mer window 的 MHC-I %rank < 2 → 風險旗標。

---

## 5. 確認事項回覆（針對 human 提的四個問題）

### Q1：沒有沈教授團隊的 PAI-1-LRP1 內部 PDB 結構

→ **使用下列公開結構**：

- **Receptor (LRP1 CR cluster)**：採用 **PDB 1J8E**（LRP1 CR cluster II，N-terminal cluster，含 CR3-CR8 tandem repeats，是文獻上 PAI-1 已知主要結合區域）。
  - 備選：2FYL（LRP1 CR5-CR6-CR7）、1CR8。
- **Full-length LRP1**：若要抓 UniProt 只為了 feature annotation，用 P98157，但只切出 CR cluster II 區段（residues 802-1184 大致範圍）做為 docking receptor。
- **PAI-1 (mimicking reference)**：PDB 3Q02（active PAI-1）當作 positive reference 序列比對，不參與 docking。
- **Peptide structure**：
  - 先丟 **AlphaFold-Multimer**（peptide + LRP1 CR cluster II）做 complex 假設。
  - 若 pLDDT_interface < 70，退而求其次用 **HPEPDOCK 2.0** blind docking 補強。
  - 最後 CHARMM/CGenFF 局部最小化降低幾何偏差（符合計畫書 p.27）。

Metadata schema 會保留 `receptor_source = "PDB:1J8E"` 以支援計畫書強調的 `evidence_source` 可追溯要求。

### Q2：有 GPU、repo 用 uv（`.venv/` 已存在）

→ **MHCflurry 2.0 安裝流程**（一次性，~500 MB 模型權重）：

```bash
uv pip install mhcflurry
uv run mhcflurry-downloads fetch models_class1_presentation
# 驗證
uv run python -c "from mhcflurry import Class1PresentationPredictor; \
    p = Class1PresentationPredictor.load(); \
    print(p.predict(peptides=['SIINFEKL'], alleles=['HLA-A*02:01']))"
```

`pyproject.toml` 目前沒有，但 `.venv/` 在，用 `uv pip install` + 之後把相依寫到 `requirements.txt`（repo 現況使用 `requirements.txt`）。

GPU 不強制需要（MHCflurry 用 TF，CPU 推論對 9-mer window sliding 已足夠快），
但 AlphaFold-Multimer / ESMFold / Rosetta FlexPepDock 會吃 GPU，這些排在 S2-S3。

### Q3：HLA panel 必要性 → 是必要的，用台灣漢人族群 panel

**為什麼必要**：MHCflurry 要指定 allele 才能預測，**全 allele 平均會稀釋風險訊號**。更重要的是，
計畫書第 8 頁明寫**「情境化校準（Contextual Calibration）…針對 NCKU PAI-1 臨床數據校準，
這是本計畫最大的數據護城河」** — 採用台灣族群 HLA panel **就是**這段論述的具體實作，
申請書審查時可以直接當佐證。

**預設 panel（台灣漢人高頻，來源：Allele Frequency Net Database + 台灣捐血中心公開資料）**：

MHC-I：
- HLA-A*11:01, A*24:02, A*33:03, A*02:07, A*02:01
- HLA-B*40:01, B*46:01, B*58:01, B*13:01, B*15:02
- HLA-C*01:02, C*07:02, C*08:01, C*03:04, C*03:02

MHC-II（給 IEDB PSSM 副指標校準參考）：
- HLA-DRB1*09:01, *15:01, *08:03, *11:01, *12:02
- HLA-DQB1*03:03, *06:01, *05:02

→ Risk aggregation：每條 9-mer window 取該 panel 內 **最強 %rank**（最壞情境），
全序列 risk = fraction of windows with `min_%rank < 2`。

### Q4：NetMHCpan 學術授權可能有 → 當內部 reference tool

策略：
- **Production / shippable**：MHCflurry 2.0（Apache-2.0，可 bundle 進 repo/SaaS）。
- **Internal reference**：若實驗室有 NetMHCpan-4.1 授權，建一個 optional adapter
  （`engine/immuno/netmhcpan_adapter.py`），僅在本地 dev 環境跑，用來**校準 MHCflurry
  的 threshold cutoff**（例如驗證 MHCflurry %rank < 2 是否對應 NetMHCpan %rank < 2 的 ~85% recall）。
- 校準結果寫成 JSON 常數落地到 repo，**不 ship NetMHCpan binary**。

---

## 6. 第一年 MVP 落地順序

| # | 步驟 | 產出 | 驗證條件 |
|---|---|---|---|
| 1 | Data Harmonization | `dataset/peptide_info/{69-80,76-88,69-88}.json` + 裸 FASTA | schema 驗證通過、保留 FITC-Ahx/MW/purity/batch metadata |
| 2 | Target Prep | `dataset/lrp1_cr2/receptor.pdb`（PDB 1J8E 切片 + CHARMM min） | CHARMM minimization 收斂、RMSD to crystal < 2 Å |
| 3 | Complex Modeling | AF-Multimer 三條 peptide × LRP1 CR2 | interface pLDDT > 70；界面殘基 ∩ 文獻 PAI-1-LRP1 接觸 > 50% |
| 4 | Docking 補強 | HPEPDOCK / ADCP 獨立 docking | Alanine mutant 預測 ΔΔG 與 `SKEMPIModel` 一致性 > 0.6 Spearman |
| 5 | Modifiable Sites List | 固定 schema JSON | schema conform、residue_index ∈ [1, N] |
| 6 | Immunogenicity | MHCflurry 2.0（台灣 HLA-I panel）+ 保留 MHC-II PSSM | 三條 peptide 的 MHC-I %rank 分布可畫；%rank<2 數量標註 |
| 7 | ADMET + Tox + Synth | ADMETlab 2.0 REST + Tox21 + AiZynthFinder | 三條 peptide 拿到完整 ADMET panel + 合成路徑評分 |
| 8 | J(x) baseline | `engine/analyzer.py` 擴充 | 初始 `w_i=0.25`；Penalty 硬條件：any window MHC-I %rank<2 or Tox21 flag |

---

## 7. 給 `engine/analyzer.py` 的具體異動清單

- `IEDBImmunoModel`（`:134-200`）：保留，作為 MHC-II 副指標。
- 新增 `MHCflurryImmunoModel`（或放到 `engine/immuno/mhcflurry_model.py`）：包裝 MHCflurry `Class1PresentationPredictor`，吃台灣 HLA-I panel，輸出 per-window %rank。
- `ProteinOptimizer.__init__`（`:326`）：同時載入兩個 immuno model。
- `score_immunogenicity`（`:374`）：改回 dict（見 §4.5）。
- 下游使用 `immunogenicity` 的三處（`:618`, `:682`, `:927`）改讀 `combined` 欄位，不影響 NSGA-II loss signature。
- `analyze` 相關輸出（`:853-866`, `:895`）：`pssm_*` 欄位保留供前端向下相容，新增 `mhc_i_*` 欄位。

---

## 8. Out of scope（第二年以後）

- **PSOS 加速演算法**（計畫書 p.7-8、p.28）：第二年於 profiling 確認瓶頸後導入。
- **子計畫 4 PBPK/QSP 回饋迴路**：需要 A2 信心分數 + 臨床暴露模擬，第三年後。
- **Causal Forest / 數位學生**：子計畫 4 的工作，不納入。
- **SRL / ReAct agent orchestration**：子計畫 2/4 的工作，子計畫 3 只需提供可被調用的 scoring API。
