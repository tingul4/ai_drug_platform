# AI-Physics 混合策略全維度藥物優化平台

**基於 SKEMPI 2.0 真實實驗數據的蛋白質藥物序列最佳化系統**

---

## 快速啟動

### macOS / Linux
```bash
chmod +x run.sh
./run.sh
```

### Windows
```
雙擊 run.bat 或在 CMD 中執行：
run.bat
```

啟動後在瀏覽器開啟 → **http://localhost:7860**

---

## 系統需求

- Python 3.9+
- 套件：`flask`, `flask-cors`, `biopython`, `numpy`, `scipy`
- 啟動腳本會自動安裝所需套件

---

## 功能說明

### 輸入
- 蛋白質氨基酸序列（單字母代碼，20-500 aa）
- 蛋白質名稱 / 靶點名稱
- 優化策略：conservative / aggressive / mixed
- 最大候選變體數（10-200）

### 輸出
- **Alanine Scan 熱點分析**：鑑定高影響力界面殘基
- **多目標 Pareto 最佳化**：同時最佳化 ΔΔG、免疫原性、聚集傾向
- **koff 動力學預測**：解離速率相對倍數
- **可藥性評分**：pI、MW、溶解度、聚集傾向
- **SKEMPI 實驗支持**：每個候選變體的實驗數據支持數量

### 資料庫
- **SKEMPI 2.0**：7,085 條真實蛋白質-蛋白質界面突變實驗記錄
- **PDB 結構**：345 個蛋白質結構（序列 + 殘基映射）
- **SQLite**：非揮發性輕量資料庫（7.2 MB），儲存分析歷史記錄

---

## 系統架構

```
ai_drug_platform/
├── db/
│   ├── skempi.db          # SQLite 資料庫（7.2 MB）
│   └── build_database.py  # 重建資料庫腳本
├── engine/
│   └── analyzer.py        # 核心分析引擎
│       ├── SKEMPIModel    # SKEMPI-calibrated 統計模型
│       └── ProteinOptimizer  # 完整最佳化流程
├── api/
│   └── server.py          # Flask REST API (port 7860)
├── frontend/
│   └── index.html         # 單頁應用前端
├── run.sh                 # Linux/macOS 啟動腳本
├── run.bat                # Windows 啟動腳本
└── README.md
```

### API 端點
| Method | Path | 說明 |
|--------|------|------|
| GET | /api/stats | 資料庫統計 |
| POST | /api/analyze | 執行最佳化分析 |
| GET | /api/history | 歷史分析記錄 |
| GET | /api/session/\<id\> | 取得特定分析結果 |
| GET | /api/skempi/search | 搜尋 SKEMPI 記錄 |
| GET | /api/pdb/\<pdb_id\> | 取得 PDB 結構資訊 |
| GET | /api/skempi/distribution | ddG 分佈（圖表用） |

---

## 重建資料庫

若需重建資料庫（例如更新 SKEMPI 數據）：

1. 將 `skempi_v2.csv` 放至 `/path/to/uploads/skempi_v2.csv`
2. 將 PDB 結構解壓至 `/path/to/skempi_pdbs/PDBs/`
3. 修改 `db/build_database.py` 中的路徑
4. 執行：`python3 db/build_database.py`

---

## 技術規格

- **ddG 預測**：SKEMPI-calibrated empirical model（340 substitution pair 統計）
- **序列比對**：k-mer Jaccard similarity（k=3）+ local alignment
- **多目標最佳化**：NSGA-II Pareto non-dominated sorting
- **koff 預測**：`koff_rel = exp(ΔΔG/RT)`，clamped to [-6, 6]
- **免疫原性**：MHC-II 9-mer sliding window scoring
- **聚集傾向**：AGGRESCAN-like hydrophobicity profiling
