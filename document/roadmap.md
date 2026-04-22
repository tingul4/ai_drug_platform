# AI 藥品優化平台 — 實作 Roadmap

依據 `document/plan.md` 與 `document/Aim3_doc_20260119.pdf` 所定義之三工具並行架構，對照 `engine/` 實作現況整理。

---

## 總覽

| 階段 | 完成度 |
|---|---|
| 架構骨架（三工具 + Pareto + API + 前端） | ✅ 已完成 |
| Tool 1：結構表徵與優化空間定義 | 🟡 僅統計層級（缺物理/AI） |
| Tool 2：動力學導向優化 | 🟡 僅 SKEMPI 查表（缺 MD/對接） |
| Tool 3：多參數成藥性/安全性評估 | 🟡 啟發式指標（缺 ADMET-AI/NetMHCpan） |
| 實驗閉環回饋機制 | ❌ 未實作 |
| 可解釋性與 LLM 報告層 | ❌ 未實作 |

---

## ✅ 已完成

### 資料與骨架層
- [x] SKEMPI 2.0 資料庫建置與清理（`db/skempi.db`, `db/build_database.py`）
- [x] IEDB 2010 MHC-II 結合資料集匯入與 PSSM 建構（`db/build_iedb_pssm.py`, `dataset/iedb_mhcii_2010/pssm.json`）
- [x] PDB 結構表與 interface residues 表
- [x] Session/Candidates 持久化（SQLite）
- [x] Flask REST API 分層（`api/server.py`）
- [x] 前端介面（`frontend/index.html`）

### Tool 1：結構表徵（統計版）
- [x] FASTA 輸入清理與驗證
- [x] k-mer + 局部比對的 PDB 相似度搜尋（`find_closest_pdb`）
- [x] SKEMPI 校準之丙胺酸掃描（`alanine_scan`）
- [x] Hotspot / warm / neutral 位點分類
- [x] BLOSUM62 近似取代分數 + 物化屬性表

### Tool 2：動力學代理（查表版）
- [x] SKEMPI 統計式 ΔΔG 預測（`SKEMPIModel.predict_ddG`）
- [x] 相對 k_off 與 residence time 估計（`predict_koff`, 由 SKEMPI koff_ratio 表）
- [x] 多策略變體產生：conservative / aggressive / mixed（`generate_variants`）
- [x] 變體多指標評分（`score_variants`）

### Tool 3：成藥性/安全性（啟發式版）
- [x] AGGRESCAN-inspired 聚集傾向（`calc_aggregation`）
- [x] CamSol-inspired 溶解度（`calc_solubility`）
- [x] pI / MW / 疏水性 / H-bond 估計
- [x] IEDB PSSM 免疫原性模型（heuristic / pssm / hybrid 三模式）
- [x] 自建 NSGA-II 非支配排序（3 目標 Pareto）

### 小分子分支（POC）
- [x] KRAS G12D POC：RDKit ETKDGv3 + Meeko + AutoDock Vina 對接
- [x] ADMET-AI（GAT）多任務預測
- [x] 描述子式 koff proxy
- [x] PyMOO NSGA-II Pareto 前緣

---

## ❌ 未完成

### Tool 1：蛋白表徵與結構預測（核心缺口）
- [ ] **ESM-2 殘基級嵌入**（Meta FAIR, L×1280）
- [ ] **AlphaFold2 / ESMFold 結構預測**（目前只做序列比對，沒有結構）
- [ ] **Uni-Mol / Uni-Mol+ 小分子 3D 表徵**（1×512 張量）
- [ ] **CHARMm 力場預處理**：氫原子補全、偏電荷分配、能量最小化
- [ ] 替代方案：OpenMM + AmberTools（開源）
- [ ] 產出**優化後 PDB 檔**（現在只回序列與指標）
- [ ] ΔG_stability（真實力場能量，而非 ΔΔG proxy）
- [ ] 局部 RMSD / local strain 指標
- [ ] 結構品質分級（實驗 vs AF2 預測）

### Tool 2：動力學真實計算
- [ ] **FoldX 物理能量估計**（突變 ΔΔG 校準）
- [ ] **ProThermDB 蛋白穩定性資料整合**
- [ ] **CDOCKER / AutoDock Vina 蛋白對接**（半柔性，蛋白端目前無對接流程）
- [ ] **AI-MD 模擬**：OpenMM + TorchMD-NET 神經網絡力場
- [ ] **OpenMMDL Flask GUI 自動化部署**
- [ ] **τ_RAMD 流程**（隨機加速力快速估 k_off 與 residence time）
- [ ] MoDEL / mdCATH MD 軌跡參考
- [ ] **3D-GNN 注意力機制**（PyTorch Geometric）
- [ ] **SHAP 歸因分析**（動力學預測可解釋性）

### Tool 3：成藥性/安全性真實工具
- [ ] **AGGRESCAN 3D**（3D 結構級聚集熱點，非目前序列啟發式）
- [ ] **CamSol** 官方工具
- [ ] **NetSurfP-3.0** 相對溶劑可及性（RSA）
- [ ] **Spatial Charge Map（SCM）** 表面電荷分佈
- [ ] **NetMHCpan-4.1**（現在以 IEDB PSSM 近似）
- [ ] **FcRn OpenMM 動態結合模擬**（抗體藥物動力學）
- [ ] **ADMET-AI / ADMETlab 2.0 蛋白端整合**（目前只接小分子）
- [ ] **Tm 預測**（熱穩定性）
- [ ] **BioPython 理化性質批次計算**
- [ ] hERG / CYP / BBB / LogP 完整小分子 ADMET 欄位對齊 PDF 規格

### 跨域決策層
- [ ] **PyMOO 正式框架**（目前自建 NSGA-II，可替換為 Platypus/PyMOO）
- [ ] **跨靶點選擇性分析**（PDBtools + RDKit 結構差異對比）
- [ ] 多目標 ≥5 維（目前只 3 維：binding / immuno / aggregation）

### I/O 規格
- [ ] 接受 **PDB 檔輸入**（目前只吃 FASTA）
- [ ] 輸出 **SDF 檔（小分子 3D 構象）**
- [ ] 輸出 **優化後 PDB 結構模型**
- [ ] Summary Table 對齊 PDF 規格（pKd / pIC50 / koff 傾向 / LogP / hERG / CYP3A4）

### 可解釋性與報告
- [ ] **3D 結構歸因圖譜**（Attention Map + SHAP 熱點疊加於結構）
- [ ] **LLM Markdown 敘事報告**（解釋 Pareto Front 化合物的藥化設計理由）
- [ ] 建議採用 **LLaMA-3 本地部署** + prompt engineering
- [ ] PyMOL / ChimeraX / py3Dmol 視覺化整合
- [ ] 數位風險標籤（Rule-based + SHAP 分級）

### 閉環驗證（Aim 3C）
- [ ] 細胞實驗資料欄位定義（IC50、細胞毒性、機制指標）
- [ ] 動物實驗資料欄位（藥效、安全性、組織分布）
- [ ] 實驗 → 計算的回饋管線（權重重分配 / 超參數調整）
- [ ] 預測 vs 實驗偏差追蹤儀表板

---

## 建議優先序

1. **Phase 1（結構化）**：ESMFold → OpenMM 能量最小化 → PDB 輸入/輸出
2. **Phase 2（動力學）**：FoldX + ProThermDB → τ_RAMD → AI-MD
3. **Phase 3（成藥性）**：AGGRESCAN 3D + NetSurfP-3.0 + NetMHCpan-4.1
4. **Phase 4（可解釋）**：3D-GNN + SHAP + LLaMA-3 報告
5. **Phase 5（閉環）**：實驗資料 schema + 回饋機制
