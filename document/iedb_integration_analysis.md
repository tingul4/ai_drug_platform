# IEDB MHC-II 2010 資料集整合分析

## 1. 整合內容摘要

| 項目 | 原始狀態 | 整合後 |
|---|---|---|
| 免疫原性模型 | MHC-II anchor 疏水熱點啟發式（`calc_immunogenicity`） | 新增 **IEDB PSSM 模型** + 保留原啟發式並提供 3 模式切換 |
| 資料依據 | 無實驗數據 | **44,541 筆 IEDB 2010 真實 MHC-II 結合親和力**（僅取人類 HLA-DR/DP/DQ，共 71,800 binder 肽 + 107,856 non-binder 肽） |
| API | 無免疫原性專用端點 | 新增 `POST /api/immuno/score`；`/api/stats` 增加 PSSM 載入狀態；`/api/analyze` 增加 `immuno_mode` 參數 |
| 前端 | 單一免疫原性欄位 | 新增「免疫原性模型」下拉（PSSM / 啟發式 / 混合），結果頁新增對照面板顯示 heur vs PSSM 與 Top-5 疑似 epitope |

## 2. 檔案變更

```
dataset/iedb_mhcii_2010/              ← 新增資料集（解壓後 5.8 MB）
├── README.txt                        ← IEDB 原始說明
├── class_II_all_split_5cv/           ← 28 alleles × 5-fold 切分
├── class_II_similarity_reduced_5cv_sep/
└── pssm.json                         ← 建立好的 9×20 PSSM 與百分位錨點

db/build_iedb_pssm.py                 ← 新增：從原始資料產生 pssm.json
engine/analyzer.py                    ← 新增 IEDBImmunoModel、calc_immunogenicity_pssm、score_immunogenicity
api/server.py                         ← 新增 /api/immuno/score、immuno_mode 透傳
frontend/index.html                   ← 新增下拉選單 + 對照面板
```

## 3. PSSM 建構方法

* **Binder 門檻**：IC50 ≤ 500 nM（MHC-II 標準閾值）
* **資料來源**：僅使用 `*_train_random_*.txt`，避免污染任何未來 benchmark 的測試切分
* **窗口處理**：每個肽（9–25 aa）被切分成所有 9-mer 子窗口，每個子窗權重 = 1/(肽長-8)，避免長肽壓倒短肽
* **分數定義**：`PSSM[pos][aa] = log2((p_binder + ε) / (p_background + ε))`
* **風險分數（0–1）**：以兩個分量混合
  * `max_component = (max_window_score - p05_binder) / (p95_binder - p05_binder)`（密度不敏感、衡量最危險 epitope 強度）
  * `density = n_hot_windows / n_windows`（衡量整條序列含 epitope 的比例）
  * `risk = 0.6 * max_component + 0.4 * density`

百分位錨點（自資料計算）：
```
binder_p05  = -0.38     binder_p50  = 1.12     binder_p95 = 2.85
nonbinder_p50 = 0.40    nonbinder_p95 = 2.16
```

## 4. 結果對照（以 Barnase 108 aa 為例）

| 模式 | 免疫原性分數 | 觀察 |
|---|---|---|
| 原始啟發式 | **0.990**（飽和） | 啟發式為「每個 9-mer 窗口加 0.01」，長序列自然累積到近上限，喪失解析度 |
| IEDB PSSM | **0.350** | 102 個窗口中僅 5 個 ≥ binder_p50；最高分 `YLQTYHKLP` (1.398)、`ITKSEAQAL` (1.345)、`HYQTFTKIR` (1.333) |
| 混合 | **0.670** | 0.5×啟發式 + 0.5×PSSM |

**結論**：啟發式在 >80 aa 蛋白幾乎恆為 0.9+，Pareto 前沿的免疫原性維度等於噪音；換成 PSSM 後，才真正能在 NSGA-II 的三目標排序中起排序作用。

## 5. 呼應 response.md 的要點

| response.md 提到的委員關切 | 本次整合如何回應 |
|---|---|
| 「胜肽路線免疫原性疑慮」（§2） | 把免疫原性從啟發式升級為 44K 筆真實 MHC-II 實驗數據建立的 PSSM，符合「候選設計與篩選階段即會納入免疫原性」的承諾 |
| 「藥動學/毒理如何映證 AI 預測精確度」（§3） | PSSM 模型可輸出具體的 Top-K epitope（含位置、序列、分數），這些是可以直接交給濕實驗端（T-cell 活化測定、HLA binding assay）做驗證的具體假設，符合「設計—模擬—回饋—再設計」中設計端產出 testable hypothesis 的角色 |
| 「整合/銜接子計畫二、四」（§4） | 免疫原性輸出為標準化欄位（risk 0–1 + epitope list + source metadata），可直接被子計畫四的 PBPK/QSP 模型讀為「免疫反應風險」輸入 |
| 子三不只跑單一親和力（原 Approach） | Pareto 多目標（ΔΔG / 免疫原性 / 聚集）現在三個維度都有實驗數據校準：ΔΔG ← SKEMPI 2.0、免疫原性 ← IEDB 2010、聚集 ← AGGRESCAN-like（未來可進一步替換） |

## 6. 後續可強化項

1. **Allele-specific 模型**：目前 PSSM 是 pan-HLA 聚合，未來可改為 allele 集合加權（依人群頻率），產出「全球/東亞/歐系」三種免疫原性風險
2. **DeImmunize 建議**：PSSM 已能計算單一氨基酸變化對 epitope 分數的影響，可延伸為「自動 epitope 消除」的變體生成策略，接入 `generate_variants` 的 `strategy` 選單
3. **ProThermDB**：dataset.md 中另一個 P0，下一輪可接入以校準 `ddG_stability`（目前只是以 0.4 權重近似 ddG_binding）
4. **ClinTox**：回應 response.md 毒性疑慮的最後一塊拼圖，可作為 P2 候補

## 7. 驗證方式

```bash
# 重建 PSSM（如果改了參數）
python db/build_iedb_pssm.py

# 啟動服務
./run.sh   # 或 python api/server.py

# 單點查詢免疫原性
curl -X POST http://localhost:7860/api/immuno/score \
  -H "Content-Type: application/json" \
  -d '{"sequence":"AQVINTFDGVADYLQTYHKL..."}'

# 完整優化（指定模式）
curl -X POST http://localhost:7860/api/analyze \
  -H "Content-Type: application/json" \
  -d '{"sequence":"...","immuno_mode":"pssm","max_variants":80}'
```

前端開啟 http://localhost:7860 後，在「優化分析」頁面可看到新的「免疫原性評分模型」下拉；執行後在「結果」頁會看到「免疫原性模型對照」面板，顯示 heuristic / PSSM 兩個分數與 Top-5 疑似 MHC-II epitope。
