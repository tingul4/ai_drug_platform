# AI-Physics Drug Optimization Platform 核心分析與預測模型

在這份系統中（以 `engine/analyzer.py` 為核心），主要使用了一個結合 **SKEMPI 2.0 實驗數據庫**的經驗模型（Empirical Model）以及多種基於物理化學屬性的啟發式（Heuristic）演算法，來進行蛋白質變異的分析與預測。

以下是具體的模型與其預測分析內容：

## 1. 核心預測模型：`SKEMPIModel`
這是一個基於真實統計數據與物理法則混合的預測模型。
* **結合自由能變化預測 (ΔΔG, `predict_ddG`)**：預測單一點突變對蛋白質結合親和力的影響。如果 SKEMPI 中有該突變的實驗數據，會直接使用統計平均值；若無，則會退回使用**基於物理特性的估算模型**（計算疏水性、電荷、體積的差異）。也會根據胺基酸位於介面、表面或內部來進行微調。
* **解離動力學預測 ($k_{off}$ ratio, `predict_koff`)**：利用 SKEMPI 的動力學數據預測突變後的相對解離速率，進而推算藥物停留時間（Residence Time）。

## 2. 成藥性與物理化學特性預測 (Drugability)
程式內建了多個經驗公式與演算法來預測蛋白質的成藥特性：
* **免疫原性風險 (`calc_immunogenicity`)**：採用類似 **NetMHCpan** 的概念，透過計算序列中是否含有 MHC-II 偏好的結合熱點（疏水性/芳香族胺基酸）來評估引發人體免疫反應的風險。
* **聚集傾向 (`calc_aggregation`)**：採用類似 **AGGRESCAN** 的算法，計算連續疏水性胺基酸的數量，評估蛋白質在溶液中產生聚集（Aggregation）的機率。
* **溶解度預測 (`calc_solubility`)**：採用類似 **CamSol** 的演算法，根據電荷、極性與疏水性比例來推算蛋白質的溶解度。
* **等電點 (pI)**：利用 Henderson-Hasselbalch 近似法推算。

## 3. 計算分析流程
* **序列相似度搜尋 (`find_closest_pdb`)**：使用快速的 K-mer Jaccard 相似度與簡化的 Smith-Waterman 對齊算法，找出 PDB 資料庫中最相似的已知結構。
* **丙氨酸掃描 (Alanine Scanning, `alanine_scan`)**：在電腦中模擬將序列上的胺基酸逐一突變為丙氨酸 (Ala)，算出 ΔΔG，以此找出對維持結合力最關鍵的「熱點 (Hot-spots)」。
* **多目標帕累托最佳化 (Pareto Optimization, `pareto_optimize`)**：使用了類似 **NSGA-II** 的非支配排序演算法 (`_fast_nds`)，對生成的變異體進行三個維度的綜合評比找出最佳解：
  1. $F_1$: 結合自由能（越低越好，代表結合力強）
  2. $F_2$: 免疫原性（越低越好，代表安全性高）
  3. $F_3$: 聚集傾向（越低越好，代表成藥性佳）

## 4. 輔助計算工具

* **分子量 (`calc_mw`)**：根據各胺基酸殘基質量加總，換算為 kDa。
* **平均疏水性 (`calc_hydrophobicity`)**：採用 Kyte-Doolittle 量表計算序列平均疏水性。
* **氫鍵估算 (`_estimate_hbonds`)**：統計序列中氫鍵供體（STNQHRKYW）與受體（STNQDEHKRY）殘基數量，以 0.6 係數估算總氫鍵數。
* **BLOSUM 替換分數 (`blosum_score`)**：根據電荷族群、極性、體積差異、疏水性差異，計算-4 ～ +4 的保守性替換分數（近似 BLOSUM62）。

## 5. 變異體生成流程 (`generate_variants`)

基於丙氨酸掃描的熱點與 warm 位置（最多取前 15 個），依三種策略生成最多 80 個變異體（另加 WT 基準）：

* **conservative**：以 BLOSUM 分數為依據，偏好保守性替換（取分數最高的 6 個候選胺基酸）。
* **aggressive**：優先挑選 SKEMPI 數據顯示 ΔΔG < −0.3 kcal/mol 的強化結合替換。
* **mixed（預設）**：60% 機率採 SKEMPI 引導，40% 採 BLOSUM 保守替換，兼顧探索與穩定性。

每個變異體的突變數以機率分布決定（40% 單突變、最多 3 個突變），隨機種子固定為 42 以確保可重現性。

## 6. 變異體評分 (`score_variants`)

對所有生成的變異體計算以下指標：

| 指標 | 說明 |
| --- | --- |
| `ddG_binding` | 各突變 ΔΔG 加總（結合自由能變化） |
| `ddG_stability` | 穩定性 ΔΔG（以 0.4 權重計） |
| `koff_relative` | 各突變 $k_{off}$ ratio 連乘（限縮至 0.001–1000） |
| `residence_time` | 藥物停留時間估算（= 10 / koff_relative） |
| 免疫原性 / 聚集 / 溶解度 / pI / MW / 氫鍵 | 物化特性完整計算 |
| `skempi_support` | 此變異體有 SKEMPI 實驗數據支撐的突變數量 |

若 PDB 介面資訊可用，會以實際介面殘基位置判斷 `is_interface`，否則退回以序列位置比例估算。

## 7. 完整優化 Pipeline (`run_optimization`)

串接上述所有步驟的主入口函式，執行順序如下：

1. 序列清洗與驗證（移除非標準胺基酸，最短 6 個殘基）
2. `find_closest_pdb` → 取得最佳 PDB 比對
3. `alanine_scan` → 找出熱點位置
4. `generate_variants` → 生成候選變異體
5. `score_variants` → 全面評分
6. `pareto_optimize` → Pareto 排序
7. 依（Pareto 排名, ΔΔG）排序，回傳前 50 名候選
8. `_save_session` → 將結果寫入 SQLite（`analysis_sessions` + `candidates` 表）

## 8. 資料持久化 (`_save_session`)

每次分析結果存入本地 SQLite（`skempi.db`），包含：

* `analysis_sessions`：session 元數據（序列、名稱、目標蛋白、執行時間、PDB 比對結果）
* `candidates`：每個候選變異體的完整指標，及一個**複合分數**（`score_composite`）：

$$\text{score} = -\Delta\Delta G_{\text{bind}} \times 0.4 + (-\ln k_{\text{off}}) \times 0.3 + (1 - \text{immunogenicity}) \times 0.15 + \text{solubility} \times 0.15$$

## 總結
此模組主要是透過**已知實驗數據驅動的統計模型**搭配**結構與物理化學演算法**，篩選並預測出親和力更強且更容易開發成實際藥物的蛋白質變異體。整體流程從序列輸入到 Pareto 最佳化結果輸出，完全在本地執行並持久化至 SQLite，無需外部 API 依賴。
