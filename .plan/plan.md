# Objective
將 AIM3 平台中胜肽分析的 Mock Data 改由實際分析工具 (Boltz-2, MMGBSA) 的 Pipeline 串接取代，以真實的工具呼叫流程搭配 Dry-run 模式驗證 3 條 PAI-1 胜肽。同時，清理遺留之舊有 SKEMPI 分析資訊。本文件為 Plan mode 生成之 ROADMAP 實作計畫。

# Key Files & Context
- `document/agent/ROADMAP.md`: 確認進度為「核心模組更新：移除舊有 SKEMPI 依賴，實作 Tool A/B 串接介面」。
- `dataset/peptide_pai1/`: 存放 3 條胜肽的基礎資料與遺留分析檔。
- `engine/analyzer.py`: 分析引擎核心。
- `api/server.py`: 長時任務 Job Queue 進入點。

# Implementation Steps

## 1. 刪除遺留分析資料 (Delete Legacy Data)
- 清除基於 SKEMPI 年初模型的 Legacy 資料，對齊 ROADMAP 中的「移除舊有 SKEMPI 依賴」：
  - 刪除 `dataset/peptide_pai1/benchmark/`
  - 刪除 `dataset/peptide_pai1/candidates/`
  - 刪除 `dataset/peptide_pai1/modifiable_sites/`
  - 刪除 `dataset/peptide_pai1/summary.json`
- 僅保留原始定義檔：`peptides.json`, `peptides.fasta`, `schema.json`。

## 2. 實作真實分析工具 Pipeline (Analyzer Updates)
- 在 `engine/analyzer.py` 實作：
  - 更新 `run_tool_a`：加入實際的 `subprocess.run(["boltz", ...])` 呼叫邏輯。
  - 更新 `run_tool_b`：加入實際的 `subprocess.run(["gmx_MMPBSA", ...])` 等指令邏輯。
  - 導入 `dry_run=True` 參數：支援快速驗證，當開啟時略過耗時計算並回傳假定或預算好的結果，以符合快速開發串接的需求。
  - 重構 `run_alanine_scan` 與 `load_candidates` 等原先依賴舊 JSON 的邏輯，改為即時產生 Dry-run 假資料或拋出尚未實作之例外以利串接。

## 3. 調整 Job Queue API (API Server)
- 於 `api/server.py` 修改：
  - 讓 `/api/jobs` 的 `create_job` 支援接收 `dry_run` 參數。
  - 將背景 Worker 與新的 `analyzer.py` 工具呼叫正確串接。

# Verification & Testing
- 測試 `/api/jobs` Endpoint，輸入 `69-80`, `76-88`, `69-88` 等三條序列。
- 帶入 `dry_run=True` 參數，驗證 Pipeline 任務狀態從 `pending` -> `processing` -> `completed`，確認 Job Queue 與狀態更新機制運作正常。