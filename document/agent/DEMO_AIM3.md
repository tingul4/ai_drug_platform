# AIM3 蛋白質藥物優化平台 — Demo 使用說明與階段性產出報告

> **⚠️ 狀態：部分過時（Partially Superseded）**
>
> 本文件記錄 2026-04-22 以前的 demo 流程，對應**通用蛋白質優化**情境
> （以 SKEMPI 查表 + IEDB-MHC-II PSSM 為核心，demo 用 hGH / TEM-1 / IFNα / Elastase / KRAS）。
>
> 在 2026-04 根據 `工程CM03-結合-unify-20260201-v35.pdf` 重新定位後，
> 子計畫三第一年主線為 **PAI-1 胜肽 → LRP1 受體** 路線，demo 與評分指標會隨之改變：
>
> - **Demo 輸入**：改以三條 PAI-1 mimicking peptides（69-80 / 76-88 / 69-88）+ LRP1 CR cluster II (PDB 1J8E) 為預設。
> - **免疫原性**：改以 MHCflurry 2.0（MHC-I，台灣 HLA panel）為主，現有 IEDB-MHC-II PSSM 降為副指標。
> - **Pareto 維度**：改為 `f_bind, f_admet, f_synth, −Penalty`（計畫書 p.28 J(x) 定義）。
>
> 本文件保留的有效部分：
> - `§1.1` 啟動方式（仍適用）
> - `§2.0` sequence / name / target 欄位關係（仍適用）
> - `§4` SQLite session 持久化（仍適用）
>
> **已過時**：`§1.2` Demo UniProt ID 清單（hGH 等）、基於 3-objective Pareto 的結果呈現。
> 現行第一年 roadmap 請見：`document/agent/roadmap.md`。

---

> 對齊文件：`document/Aim3_doc_20260119.pdf`
> 報告日期：2026-04-22

---

## 1. Demo 快速上手

### 1.1 啟動

```bash
cd /raid/danielchen/ai_drug_platform
./run.sh                     # Linux / macOS
# 瀏覽器開啟 http://localhost:7860
```

### 1.2 Demo 輸入目標（透過 UniProt ID 一鍵載入）

| UniProt ID | 蛋白名稱 | 為什麼選它 |
|---|---|---|
| `P01241` | Human Growth Hormone (hGH) | SKEMPI 含 251 條 hGH / hGHbp 界面突變資料，Pareto 優化結果豐富 |
| `P62593` | TEM-1 β-lactamase | SKEMPI 277 條 vs BLIP，抗生素抗藥性研究經典標的 |
| `P01563` | Interferon α-2 | SKEMPI 220 條，cytokine-receptor 介面突變完整 |
| `P08246` | Human Leukocyte Elastase | SKEMPI 259 條，serine-protease / inhibitor 典型案例 |
| `P01116` | KRAS | 對應小分子 POC (PDB 7RPZ)，可搭配小分子管線 demo |

呼叫方式（後端 API）：

```bash
curl http://localhost:7860/api/uniprot/P01241        # 取序列 + 特徵 + PDB/AF 交叉參考
curl http://localhost:7860/api/uniprot/demo          # 列出以上 demo 清單
```

回傳範例（節錄）：
```json
{
  "id": "P01241",
  "protein_name": "Somatotropin",
  "length": 217,
  "sequence": "MATGSRTSLL...",
  "features": [
    {"type": "Signal",      "start": 1,   "end": 26},
    {"type": "Domain",      "start": 33,  "end": 217, "description": "Somatotropin"},
    {"type": "Binding site","start": 64,  "end": 64}
  ],
  "pdb_xrefs": ["1A22", "1AXI", "3HHR", ...],
  "alphafold_id": "AF-P01241-F1"
}
```

### 1.3 執行分析

打 `/api/analyze`（或用前端表單），核心 payload：

```json
{
  "sequence":     "FPTIPLSRLFDNAMLRAHRLHQLAFDTYQ...",
  "name":         "Somatotropin",
  "target":       "hGH receptor",
  "strategy":     "mixed",
  "max_variants": 80
}
```

目前 `/api/analyze` 回傳的是**突變序列集合與多目標評分**（未輸出 3D 結構檔，照本階段計畫）。

---

## 2. 參數選項與意義

### 2.0 三個輸入欄位到底是什麼關係？（常見誤解）

平台目前的輸入可以想成 **「一個藥（會被突變的蛋白）」 + 「它要打的標靶」** 兩個物件：

| 欄位 | 角色 | 目前實作 | 備註 |
|---|---|---|---|
| `sequence` | **藥物蛋白的氨基酸序列** — 就是被優化、會被打突變的那個分子 | ✅ 真正進入 ΔΔG 模型運算 | 例：抗體 Fab、胜肽藥、工程化細胞因子 |
| `name` | 藥物名稱 | ⚠️ 只當 session 標籤，不影響運算 | 方便在歷史記錄裡辨識 |
| `target` | **標靶蛋白名稱** — 例：hGH receptor、IL-2Rα | ⚠️ 目前只是字串標籤 | 標靶資訊是透過 SKEMPI 找到的 PDB 複合物「間接」帶入（每條 SKEMPI 記錄本來就是蛋白 A–蛋白 B 的界面） |

**對應 PDF 的概念：**
- PDF「候選藥物數據」 = 這裡的 `sequence`（+ 未來要接的藥物 PDB）
- PDF「標靶與生物環境數據」 = 這裡的 `target`（下一階段會改成 `target_uniprot_id` + 拉標靶 PDB/AF 結構）

> **下一階段升級建議的 API 形狀：**
> ```json
> {
>   "drug_uniprot_id":   "P01241",    // 或 drug_sequence
>   "target_uniprot_id": "P10912",    // 或 target_sequence / target_pdb_id
>   "strategy": "mixed", "max_variants": 80, "immuno_mode": "pssm"
> }
> ```

### 2.1 分析選項

| 參數 | 選項 | 含義與建議 |
|---|---|---|
| `sequence` | FASTA 字串 | 藥物蛋白序列；6–2000 aa（超出會回 error） |
| `name` | 自由字串 | session 標籤 |
| `target` | 自由字串 | 標靶名稱（目前只是標籤） |
| `strategy` | `conservative` / `aggressive` / `mixed` | 突變策略（見 §2.2） |
| `max_variants` | 10–200 | 候選變體數量；80 為 demo 合理值 |

> 免疫原性一律以 IEDB-2010 PSSM 計算，目前不提供 `immuno_mode` 參數（heuristic / hybrid 為下階段選項；見 §7）。

### 2.2 突變策略（strategy）

- **`conservative`**：BLOSUM62 相近性替換（同電荷、同極性、體積差小）。安全；維持折疊；ΔΔG 變動小；適合提親和力微調或做 humanization。
- **`aggressive`**：從 SKEMPI 統計中挑「此位替換後平均 ΔΔG_binding < −0.3 kcal/mol」的氨基酸，傾向大幅強化結合；風險：可能破壞穩定性或引入免疫原性。
- **`mixed`**（預設）：60 % SKEMPI-guided（資料驅動）+ 40 % 保守（BLOSUM 保底）。是目前提交報告建議的預設。

### 2.3 免疫原性計分（IEDB-2010 PSSM）

目前僅實作 **PSSM** 一種模式，由 `dataset/iedb_mhcii_2010/class_II_all_split_5cv`（IEDB 官方 MHC-II 2010 benchmark）建出 9-mer PSSM：
- 訓練規模：**71,800 個 binder peptides / 107,856 個 non-binder peptides**（由原始 44,541 筆 binding-affinity 測量、28 個 HLA-DR/DP/DQ alleles 經 split + ≥ 9-mer 展開得到）
- 計分方式：每個 9-mer 窗口算 log-odds，取「max-window 分數 × 0.6 + 超過 binder 中位數的窗口比例 × 0.4」→ `risk ∈ [0, 1]`
- 產出欄位：`immunogenicity`、`pssm_top_hits`（前 5 個高風險片段）、`pssm_n_hot`、`pssm_n_windows`

> 下階段會補 `heuristic`（MHC-II 錨點疏水/芳香殘基計數）與 `hybrid`（兩者 0.5/0.5 混合）做 sanity-check；見 §7。

---

## 3. 為什麼用「丙胺酸（Alanine）掃描」？

**Alanine scanning** 是界面熱點識別的黃金標準方法（Clackson & Wells, 1995）：

1. **Alanine 的側鏈最小（只有 −CH₃）但不破壞主鏈**。把某位殘基換成 A 時，等於「移除側鏈但保留骨架」，讓我們量化該殘基側鏈對結合能的貢獻。
2. **Δ(結合自由能) = ΔG(mutant) − ΔG(WT)** 若 ≥ 1 kcal/mol，該殘基被定義為**熱點（hot-spot）**。這是 SKEMPI 資料庫的標準閾值。
3. **熱點通常佔界面 10–20 % 殘基，卻貢獻 > 80 % 的結合能**。聚焦突變這些位點，搜尋空間從 `20^L` 壓到 `20^(0.15L)`，計算量下降數個數量級。

本平台 `engine/analyzer.py::alanine_scan()` 對序列每個位點預測 ΔΔG(X→A)，用 SKEMPI 統計做 calibration，並標註 `is_hotspot`（閾值 1.0 kcal/mol）。後續突變庫只在熱點+warm 位點（top-15）採樣。

---

## 4. 主要觀察指標（Output Metrics）

### 4.1 候選變體 Summary Table — 對齊 PDF Page 4

PDF 的 Summary Table 原始設計為小分子（`ID | pIC50 | pKd | koff | LogP | hERG | CYP3A4`）。
蛋白質藥物沒有 LogP / hERG / CYP3A4，但 PDF Page 4 蛋白質段落指定了對應的「去風險」欄位。
本平台的候選表頭依此對齊：

| PDF 欄位 | 小分子意義 | **本平台（蛋白質）對應欄位** | 計算方式 |
|---|---|---|---|
| ID | 候選編號 | `variant_id`（V0001…） | 生成時指派 |
| pIC50 (活性) | 活性預測 | `pIC50_abs`（有 PDBbind 錨點時）／ `pIC50_shift`（無錨點） | WT pKd + `−ΔΔG/1.363` |
| pKd (親和力) | 親和力 | `pKd_abs`（有 PDBbind 錨點時）／ `pKd_shift`（無錨點） | WT pKd + `−ΔΔG/(RT·ln10)` @ 298 K |
| koff 傾向 | 動力學 | `koff_relative` | `exp(ΔΔG/RT)`，< 1 = 解離更慢 |
| LogP (理化) | 小分子理化 | `solubility` | CamSol-inspired，0–1 |
| hERG (心毒性) | 小分子安全 | `immunogenicity` (MHC-II) | IEDB-2010 PSSM，0–1 |
| CYP3A4 (代謝) | 小分子代謝風險 | `aggregation` | AGGRESCAN-like，0–1 |

> **絕對 pKd / pIC50 已啟用**（PDBbind v2020 R1 整合完成）：
> 當匹配到的 SKEMPI PDB 在 PDBbind 有實驗 Kd/Ki/IC50 時（345 個 SKEMPI PDB 中有 **229 個** 覆蓋到），
> 平台會回傳 `pKd_abs = pKd_WT + ΔpKd`。無錨點時退回顯示 `Δ shift`。
> 範例（mature hGH → 1A22）：PDBbind 記錄 Kd = 0.34 nM → pKd_WT = 9.47；最佳候選 ΔpKd = +1.20 → **pKd_abs = 10.67**（≈ 21 pM，16× 強於 WT）。
>
> **PDBbind tie-breaker（2026-04-22 加入）：** `find_closest_pdb()` 若 top-1 PDB 無 PDBbind 錨點，會在 similarity 差距 ≤ 0.15 內自動 promote 有錨點的候選。例如 mature hGH 原本 top-1 是 1BP3 (sim 0.87, 無覆蓋)、top-2 是 1A22 (sim 0.79, 有覆蓋)；tie-breaker 後 1A22 升為 top-1，`pKd_abs` 成功產出。容忍度由 `ProteinOptimizer.PDBBIND_TIEBREAK_TOL` 控制。

### 4.2 完整輸出欄位

對齊 PDF Page 4 「預期輸出」：

| 指標 | 實作狀態 | 含義 | 讀值方向 |
|---|---|---|---|
| `ddG_binding` (ΔΔG) | ✅ | 與標靶結合自由能變化（SKEMPI-calibrated；多突變累加會 clamp 至 ±10 kcal/mol，見 §4.3） | 越負 = 結合越強 |
| `pKd_shift` (ΔpKd) | ✅ | 由 ΔΔG 推導的親和力變化（正 = 更強） | 越正越好 |
| `pIC50_shift` (ΔpIC50) | ✅ | 競爭性抑制時 ≈ ΔpKd | 越正越好 |
| `pKd_abs` | ✅ 新增 (PDBbind) | 絕對 pKd = pKd_WT + ΔpKd（有 229/345 SKEMPI PDB 覆蓋） | 越大越強 |
| `pIC50_abs` | ✅ 新增 (PDBbind) | 絕對 pIC50 | 越大越強 |
| `ddG_stability` | ⚠️ 粗估 | 蛋白單體摺疊穩定度；下階段接 FoldX / ProThermDB；**已納入 Pareto (F4)** | 越負越穩定 |
| `koff_relative` | ✅ | 解離速率相對值 `exp(ΔΔG/RT)` | < 1 = 解離更慢（好） |
| `residence_time` | ✅ | `10 / koff_rel`，動力學停留時間比例 | 越大越好 |
| `immunogenicity` | ✅ (IEDB PSSM) | T-cell 表位風險 0–1 | 越低越安全 |
| `aggregation` | ✅ (AGGRESCAN-like) | 聚集傾向 0–1 | < 0.4 可接受 |
| `solubility` | ✅ (CamSol-inspired) | 溶解度 0–1 | 越高越好 |
| `pI` | ✅ | 等電點，影響配方 pH | 建議離生理 7.4 至少 ±1 單位 |
| `mw_kDa` | ✅ | 分子量，做 PK 預估 | — |
| `skempi_support` | ✅ | 該突變在 SKEMPI 有多少條實驗證據 | 越多越可信 |
| `pareto_rank` | ✅ | 多目標非支配排序層級（**5 目標**：ddG_binding / immuno / agg / ddG_stability / koff） | 1 = 最優前緣 |
| `warnings` | ✅ 新增 (§4.3) | 物理值範圍警告標籤（ddG_extrapolated / pkd_capped / aggregation_high…） | 陣列越短越好 |
| **`Tm`（熱穩定）** | ❌ 未實作 | 見 §7（優先序 #3），由 `ddG_stability` 粗估頂著 | — |
| **優化後 PDB 3D** | ❌ 本階段略 | 本階段只回傳突變序列 | — |
| **Attention / SHAP 熱圖** | ❌ 本階段略 | PDF Page 5 可解釋性素材 | — |
| **LLM 敘事報告** | ❌ 本階段略 | PDF Page 5 執行摘要 | — |

### 多目標 Pareto 優化（2026-04-22 擴成 5 目標）

同時最小化五個目標：

```text
F1 = ddG_binding        (親和力 — 越負越好)
F2 = immunogenicity     (免疫原性 — 越低越好)
F3 = aggregation        (聚集傾向 — 越低越好)
F4 = ddG_stability      (單體摺疊穩定性 — 越負越好；目前為粗估，Tm 接上後替換)
F5 = koff_relative      (相對解離速率 — 越小越慢 = 停留時間越長)
```

使用 NSGA-II 快速非支配排序產生 Pareto front。`pareto_rank=1` 的候選即為「沒有其他候選能在這五個目標上同時更好」的最優集合，提交給化學家做下一步實驗驗證。

> **為什麼沒把 pKd 放進 Pareto？**
> `pKd_shift = −ΔΔG_binding / 1.363` 是 ΔΔG 的線性變換，相關係數 ≈ −1.0，納入等於把親和力軸乘兩次、稀釋其他目標。pKd / pIC50 仍以展示欄位呈現（見 §4.2）。

> **為什麼沒把 solubility / pI 放進 Pareto？**
> 這兩者是「通過/不通過」的配方學閾值（見 §4.3），放進 Pareto 只會塞爆 front；改以警告標籤（`warnings`）呈現，候選表格可直接篩掉。

---

### 4.3 物理值範圍與警告閾值（2026-04-22 新增）

為避免模型外推或候選超出物理上合理範圍，每個候選（含 WT baseline）都會帶一個 `warnings: [...]` 欄位，前端/JSON 可直接讀。clamp 會**修改數值**；warning 只**加標籤不改值**。

| 欄位 | 合理範圍 | 超出時動作 | 標籤 |
|---|---|---|---|
| `ddG_binding`（多突變累加） | [−10, +10] kcal/mol | **clamp** 到邊界 | `ddG_extrapolated` |
| `pKd_abs` | ≤ 14（≈ fM affinity ceiling） | **cap** 到 14 | `pkd_capped` |
| `koff_relative`（單步） | [0.001, 100] | clamp（舊有） | — |
| `koff_relative`（多突變累乘） | [0.001, 1000] | clamp（舊有） | — |
| `abs(pI − 7.4)` | ≥ 1.0 pH 單位 | 僅加警告 | `pI_near_physiological` |
| `aggregation` | < 0.4 | 僅加警告 | `aggregation_high` |
| `solubility` | > 0.4 | 僅加警告 | `solubility_low` |
| `mw_kDa` | 0.5 – 200 kDa | 僅加警告 | `mw_out_of_range` |
| `len(sequence)` | 6 – 2000 aa | **拒絕** 並回傳 error | — |

**回傳範例（節錄）：**
```json
{
  "variant_id": "V0023",
  "ddG_binding": -10.0,
  "pKd_abs":     14.0,
  "aggregation": 0.52,
  "warnings":   ["ddG_extrapolated", "pkd_capped", "aggregation_high"]
}
```

> 實作位置：`engine/analyzer.py::ProteinOptimizer._physical_warnings()`，類常數（`DDG_BIND_CAP`, `PKD_ABS_CAP`, `AGG_MAX_OK`, `SOL_MIN_OK`, `MW_MIN_KDA`/`MAX_KDA`, `PI_MIN_DIST`）可集中調整。

---

### 4.4 Demo 實跑結果（mature hGH, 2026-04-22 基準）

輸入：UniProt **P01241**（mature hGH，signal peptide 切除後 191 aa，即 `sequence[26:]`），預設參數 `strategy=mixed`、`max_variants=80`。

| 指標 | 數值 |
|---|---|
| runtime | **0.24 秒**（CPU，single process） |
| matched PDB | **1A22** (hGH–receptor complex, sim 0.79)；經 PDBbind tie-breaker 由 1BP3 (sim 0.87, 無覆蓋) 降為次佳 |
| WT 實驗 Kd (PDBbind) | 0.34 nM → **pKd_WT = 9.47** |
| 熱點 (hotspots) | **22 個**（ΔΔG(X→A) > 1 kcal/mol 且位於界面） |
| 候選數 | 80 生成 / 保留 top 50 |
| Pareto rank-1 | **13 個**（5 目標下的非支配最優集合） |
| 最佳候選 | `V0023` — 雙突變 `E56H + H18T` |
| 最佳 ΔΔG_binding | **−2.10 kcal/mol** |
| 最佳 ΔpKd | **+1.54** |
| **最佳 pKd_abs** | **11.01**（≈ 9.8 pM，比 WT 0.34 nM 強 **35×**）|
| WT 警告 | `pI_near_physiological`, `aggregation_high`, `solubility_low` |
| 候選警告分佈 | `pI_near_physiological` 50/50；`aggregation_high` 50/50 |

**可講的重點：**
1. **小於 1 秒出結果**：整套管線（含 5 目標 NSGA-II、80 候選評分）在 191 aa 輸入下 0.24 秒完成；適合互動式 demo。
2. **絕對親和力可交付**：不只給 ΔΔG，而是給 **pKd_abs = 11.01**（pM 級），化學家可以直接拿這個數字對標競品。
3. **熱點集中**：191 aa 中只有 22 個位點貢獻 > 1 kcal/mol，佔 11.5%（符合 Clackson-Wells「熱點佔界面 10–20%」的經典觀察）。
4. **Pareto 壓縮效果**：80 個候選中有 13 個 rank-1，代表 5 目標同時最優的解空間約 26%，留給化學家足夠選擇但不爆炸。
5. **誠實顯示系統限制**：WT 本身就 `aggregation_high`（hGH 確實是高疏水性蛋白）；警告標籤讓評審一眼看出「這是底子就這樣，不是 false-negative」。

> 重現：`curl -X POST localhost:7860/api/analyze -H 'Content-Type: application/json' -d '{"sequence": "<mature hGH 191 aa>", "strategy": "mixed", "max_variants": 80}'`

---

### 4.5 邊界情況與已知限制

| 情境 | 系統行為 |
|---|---|
| 輸入序列在 SKEMPI 345 個 PDB 裡找不到 > 0.05 k-mer 相似度的結構 | `pdb_matches = []`、退回純 physics-based ΔΔG（`_physics_estimate`）；`skempi_support=0`、`confidence=0.2` |
| 輸入序列 < 6 aa 或 > 2000 aa | 回傳 error，不執行管線 |
| 匹配到的 PDB 不在 PDBbind（見 §4.1 hGH 範例） | `pKd_abs=None`，只顯示 `pKd_shift`；不影響 Pareto |
| 多突變 ΔΔG 累加超過 ±10 kcal/mol | clamp 到邊界，候選標 `ddG_extrapolated` 警告 |
| `(wt_aa, mut_aa)` 在 SKEMPI 統計找不到 | 退回 physics-based 估計 + `n_support=0` + 低置信度 |
| 非標準 AA（B/J/O/U/X/Z） | 清洗階段直接過濾，不會進入計算 |

---

## 5. 目前已完成的功能清單

### 5.1 核心管線（`engine/analyzer.py`）
- [x] FASTA 序列清洗與驗證（6–2000 aa）
- [x] k-mer + local-align 相似度搜尋最接近的 SKEMPI PDB
- [x] Alanine scan 熱點識別（SKEMPI-calibrated ΔΔG）
- [x] 三種突變庫產生策略（conservative / aggressive / mixed）
- [x] 多指標打分（ΔΔG、ΔpKd、pKd_abs、koff、免疫原性、聚集、溶解度、pI、MW）
- [x] **NSGA-II Pareto 非支配排序（5 目標：ddG_bind / immuno / agg / ddG_stability / koff）**
- [x] IEDB-PSSM 免疫原性模型（28 alleles → 71,800 binder peptides）
- [x] **PDBbind v2020 R1 整合**：pKd_abs / pIC50_abs（229/345 SKEMPI PDB 有覆蓋）
- [x] **物理值範圍警告（`warnings`）**：8 種標籤 + ddG/pKd clamp（見 §4.3）
- [x] Session / candidates 寫入 SQLite 歷史
- [x] 綜合分數 `score_composite` = `-0.4·ΔΔG + 0.3·(-log koff_rel) + 0.15·(1-immuno) + 0.15·solubility`

### 5.2 API 層（`api/server.py`）
- [x] `GET  /api/stats` — 資料庫統計
- [x] `POST /api/analyze` — 主分析管線
- [x] `GET  /api/history?limit=N` — 最近 N 個 session 列表
- [x] `GET  /api/session/<id>` — 讀取歷史 session
- [x] `GET  /api/skempi/search` — SKEMPI 原始記錄查詢
- [x] `GET  /api/pdb/<pdb_id>` — PDB 結構資訊
- [x] `GET  /api/skempi/distribution` — ΔΔG 全庫分佈（圖表用）
- [x] `POST /api/immuno/score` — 單序列免疫原性快速評估
- [x] **`GET  /api/uniprot/<acc>`** — UniProt 序列 + 特徵 + PDB/AF 交叉參考
- [x] **`GET  /api/uniprot/demo`** — 列出 demo 標的
- [x] `GET  /api/smallmol/poc` — 小分子 POC 結果（KRAS G12D / PDB 7RPZ）
- [x] `GET  /api/smallmol/plot` — 小分子 POC 視覺化 PNG

### 5.3 小分子 POC（`dataset/kras_g12d_poc/`）
- [x] RDKit ETKDGv3 3D 建構
- [x] Meeko PDBQT 準備
- [x] AutoDock Vina 1.2.3 對接
- [x] ADMET-AI（Graph Attention Network）ADMET 預測
- [x] 描述子-based koff proxy
- [x] PyMOO NSGA-II 非支配排序

### 5.4 前端（`frontend/index.html`）
- [x] 序列輸入、策略/參數選擇
- [x] Pareto 前緣、Alanine scan 熱圖、候選列表
- [x] 小分子 POC 展示
- [ ] UniProt ID 下拉（本次新增 endpoint，UI wiring 列 §7 待辦）

---

## 6. 使用的資料庫／工具 — 為什麼用它們？

對齊 PDF Page 7–13「可用資料庫及預測工具」。

### 6.0 常見誤解：SKEMPI 能不能用 UniProt 取代？

**不能。兩者角色不同，是互補關係：**

|  | SKEMPI 2.0 | UniProt |
|---|---|---|
| 性質 | 實驗熱力學資料集 | 序列 / 功能註解資料庫 |
| 核心內容 | 7,085 條 PPI 突變 + 實驗 ΔΔG + koff + 溫度 + 界面分類 | > 2.5 億條序列 + Domain / Active site / PTM / 疾病變異 / PDB·AF 交叉參考 |
| 在平台的角色 | **校正 ΔΔG 預測模型**（沒它模型退回純物理化學估計） | **查詢入口**：UniProt ID → 自動拉序列 + 結構線索 |
| 輸入 / 輸出關係 | 提供「突變了會怎樣」的訓練答案 | 提供「要優化哪個蛋白」的資訊 |

結論：**保留 SKEMPI 做模型校正，加上 UniProt 做輸入查詢**（本階段已完成）。

### 6.1 實驗驗證資料庫（Experimental — 高可信度）

| 資料庫 | 用途 | 規模 / 內容 | 為什麼是它 |
|---|---|---|---|
| **SKEMPI 2.0** | ΔΔG / koff 模型校正 | 7,085 條 PPI 界面突變；每條含 Kd(WT/Mut)、koff、溫度、界面分類 (COR/RIM/SUP) | PPI 突變熱力學的**標準訓練集**（CC BY 4.0）；PDF 明確列入 |
| **PDB** | 結構 / 界面註解 | X-ray、NMR、Cryo-EM 實驗結構 | 結構生物學權威來源；本平台取其界面殘基 |
| **IEDB 2010 MHC-II** | 免疫原性 PSSM | 原始 44,541 筆 binding affinity / 28 HLA-DR/DP/DQ alleles → 展開為 **71,800 binder / 107,856 non-binder** 9-mer peptides 訓練 | 免疫表位的黃金標準；CC BY 4.0 |
| **UniProt**（前次接入） | 標靶序列與功能註解 | > 2.5 億條序列；含 Domain / Active site / Binding site / PTM / 疾病變異 | PDF 明確推薦；單一 ID 即可拉回 canonical sequence + PDB / AF 交叉參考；後續可用於**保護 active-site 不被亂突變** |
| **PDBbind v2020 R1**（本次接入） | 絕對 pKd / pIC50 錨點 | 19,037 PL + 2,798 PP 複合物實驗 Kd/Ki/IC50 | 解鎖 PDF Summary Table 的絕對 pKd 欄；229/345 SKEMPI PDB 有覆蓋 |

### 6.2 預測 / 生成資料庫

| 資料庫／工具 | 用途 | 為什麼 |
|---|---|---|
| **AlphaFold-DB** | 備援結構（下一階段） | 當 SKEMPI 內 345 個 PDB 無相似匹配時，用 AF2 預測結構補位；CC BY 4.0 |
| **ADMETlab 2.0 / ADMET-AI** | 小分子 ADMET | 80+ 指標（LogP、溶解度、hERG、CYP450、BBB）；已整合在小分子 POC |

### 6.3 目前**尚未接入**但 PDF 建議（§7 待辦）

| 資料庫 | 將來補什麼 |
|---|---|
| **ProThermDB** | Tm 與單體 ΔΔG_stability 模型訓練資料（目前 stability 只是 binding 模型套殼） |
| **BindingDB** | 擴充 PDBbind 沒覆蓋的 PDB 當 WT 錨點（再補約 116 個 SKEMPI PDB） |
| **ChEMBL / DrugBank** | 小分子 SAR 與類藥分子參考 |
| **NetMHCpan-4.1** | PDF 指定的黃金標準免疫原性工具，取代現行 IEDB-PSSM |
| **AGGRESCAN 3D / Zyggregator** | 結構-level 聚集預測（目前是序列-level） |
| **FoldX / Rosetta** | 單體穩定性 ΔΔG 精確計算 |

### 6.4 資料授權一覽

| 資料來源 | 授權 | 商業使用 |
|---|---|---|
| SKEMPI 2.0 | CC BY 4.0 | ✅ 允許（須標註來源） |
| IEDB 2010 MHC-II | CC BY 4.0 | ✅ 允許 |
| UniProt | CC BY 4.0 | ✅ 允許 |
| PDB 結構檔 | 公眾領域（RCSB 條款） | ✅ 允許 |
| **PDBbind v2020 R1** | **學術免費**；商業須向中科院物理所申請授權 | ⚠️ 商業受限 |
| AlphaFold-DB | CC BY 4.0 | ✅ 允許 |
| ADMETlab 2.0 / ADMET-AI | 學術免費 | ⚠️ 商業須確認 |

> 對外 demo 若涉及商業洽談，需移除或替換 PDBbind 的絕對 pKd 功能，或取得 PDBbind 商業授權。

---

## 7. 下一階段急需補上（優先序）

依 PDF 預期輸出對齊後，**最需要補的**依序：

1. **AlphaFold-DB 備援結構管線**
   - 為什麼急：目前 SKEMPI 只有 345 PDB，UniProt ID 被接入後會有大量標的沒有對應的實驗結構。
   - 工作量：小（用 UniProt → AF ID → 抓 PDB 檔即可）。

2. **標靶蛋白輸入從字串升級為 UniProt ID + 結構**
   - 為什麼急：目前 `target` 只是字串，沒有進入計算；接入後才有「跨標靶選擇性分析」的基礎。
   - 工作量：中（需修改 `run_optimization` 讀目標結構、做界面判定）。

3. **Tm 熱穩定性預測**（本次明確列為待辦；Demo 前不實作）
   - 為什麼急：PDF 蛋白質側「去風險三面向」其一，目前是唯一的 gap（親和力 / 免疫原性 / 聚集 / 穩定性四者，穩定性靠 `ddG_stability` 粗估頂著）。
   - 目前替代：`ddG_stability` 已納入 Pareto F4（見 §4 多目標段），以 binding 模型套殼近似 — 接上 Tm 後直接替換 F4。
   - 方案選項：(a) heuristic ~1 hr（AA composition + 二硫鍵 + 疏水電荷經驗式，精度 ±8°C），(b) FoldX 單體能量 0.5–1 天，(c) ESM-2 + ProThermDB fine-tune 3–5 天。

4. **絕對親和力 (pKd / pIC50)**
   - 方案：PDBBind 訓練一個 GNN 回歸頭，把相對 ΔΔG 轉為絕對 −log(Kd)。

5. **把 Alanine scan 與 UniProt feature track 對齊**
   - Active site / Binding site 不應被當作候選突變位點（避免設計出 loss-of-function 或已知致病變異）。
   - 本次已把 features 拉進來，但 `alanine_scan()` 還沒消費這個欄位。

6. **前端 UniProt ID 下拉**
   - 小工作；把 `/api/uniprot/demo` 接到 UI 即可。

---

## 8. 一分鐘 Demo 腳本（適合口頭報告）

> 「這個平台輸入一條蛋白質藥物序列，輸出一組經過多目標優化、可進實驗的突變候選。
>
> 我們用 SKEMPI 2.0 的 7,085 筆真實 PPI 突變訓練 ΔΔG 模型；用 IEDB 2010 的 4.4 萬筆 MHC-II binding（展開 71.8k binder peptides）訓練免疫原性 PSSM；用 PDBbind v2020 的 21,835 個複合物提供絕對 pKd 錨點。
>
> 流程是：UniProt 拉序列 → Alanine scan 找熱點 → 三策略產生突變庫 → 打 9 個指標分數 → **NSGA-II 在 5 個目標上做 Pareto 排序**（親和力 / 免疫原性 / 聚集 / 穩定性 / 動力學）+ 自動附加物理值範圍警告。
>
> PDF 預期輸出對齊狀況：ΔΔG ✅、絕對 pKd ✅（229/345 覆蓋）、koff ✅、免疫原性 ✅、聚集 ✅、溶解度 ✅、Pareto ✅；還缺 **Tm、優化後 3D 結構、Attention/SHAP 熱圖、LLM 敘事報告**，排進下階段（§7）。」

---

## 9. 如何從零重現（Build Order）

本 repo 不含重建後的大型 DB（`db/skempi.db` 約 數十 MB），若要從原始資料重跑，依序執行：

```bash
# 0. 環境
pip install -r requirements.txt

# 1. 建 SKEMPI + PDB 結構主表（產出 db/skempi.db：skempi_mutations, pdb_structures, pdb_residues）
python db/build_database.py

# 2. 建 PDBbind 絕對親和力索引（產出 db/skempi.db::pdbbind_affinity，21,835 條）
python db/build_pdbbind.py

# 3. 建 IEDB-MHCII PSSM（產出 dataset/iedb_mhcii_2010/pssm.json）
python db/build_iedb_pssm.py

# 4. 啟動 API + 前端
./run.sh                         # → http://localhost:7860
```

**原始資料需自行下載並放到以下路徑**（授權見 §6.4）：

| 來源 | 放置位置 |
|---|---|
| SKEMPI v2.csv + SKEMPI_v2 PDB 檔 | `dataset/skempi_v2/` |
| IEDB 2010 MHC-II (`class_II_all_split_5cv/`) | `dataset/iedb_mhcii_2010/` |
| PDBbind v2020 R1 index | `dataset/PDBbind_v2020_R1/index/` |

> 整個 build 流程在一般 workstation 約 5–10 分鐘，其中 SKEMPI PDB 解析最耗時。Build 完 `db/skempi.db` 為 **單一 SQLite 檔**，可直接打包移動。

---

## 附錄 A：DB Schema（精要）

| 表 | 主要欄位 | 寫入時機 |
|---|---|---|
| `skempi_mutations` | pdb_id, mutation_clean, ddG_kcal, koff_ratio, location (COR/RIM/SUP) | build 時一次寫入 |
| `pdb_structures` | pdb_id, sequences (JSON per chain), proteins, resolution | build 時 |
| `pdb_residues` | pdb_id, chain, seq_num, is_interface | build 時 |
| `pdbbind_affinity` | pdb_id, complex_type (PP/PL), kind (Kd/Ki/IC50), pKd | build 時 |
| `analysis_sessions` | session_id, input_sequence, matched_pdb, runtime_sec, result_summary | 每次 `/api/analyze` |
| `candidates` | session_id, variant_id, mutations, ddG_binding, pKd_abs, pareto_rank, score_composite ... | 每次 `/api/analyze`（top 100） |

> ⚠️ **`warnings` 欄位目前只在 JSON 回傳，不寫入 `candidates` 表**。若需歷史可查詢警告，下階段需加欄位 migration。
