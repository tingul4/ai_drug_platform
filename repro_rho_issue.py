
import numpy as np
from scipy.stats import spearmanr

def simulate_rho_analysis():
    # 實驗基準 (Stefansson 1998): 重要性 K69 > K80 > K88
    # 假設能量貢獻 (越高代表突變後損失越多，即越重要)
    ground_truth = {"K69": 10.0, "K80": 5.0, "K88": 2.0}
    
    peptides = {
        "PE1016 (WT-like)": ["K69", "K80", "K88"],
        "PE1018":           ["K69", "K80"],
        "PE1020 (Missing K69)": ["K80", "K88"]
    }

    print("=== 實驗 A: 純序列模型 (不考慮結構偏移) ===")
    for name, residues in peptides.items():
        # 模擬預測：始終按照殘基本身的性質給分
        preds = [ground_truth[r] for r in residues]
        actual = [ground_truth[r] for r in residues]
        rho, _ = spearmanr(actual, preds)
        print(f"{name}: Predicted Rank={[r for r in residues]}, rho={rho:.3f}")

    print("\n=== 實驗 B: 結構敏感模型 (模擬 PDF 邏輯) ===")
    for name, residues in peptides.items():
        if "K69" not in residues:
            # 模擬 PE1020 的情況：缺少錨點導致 Pose 錯誤，排名反轉
            # 原本應是 K80(5.0) > K88(2.0)，因 Pose 錯誤變成 K88 > K80
            preds = [1.0, 4.0] # 錯誤的排名
            actual = [5.0, 2.0] # 實驗基準
            note = "(Non-native Pose detected!)"
        else:
            preds = [ground_truth[r] for r in residues]
            actual = [ground_truth[r] for r in residues]
            note = "(Native Pose)"
        
        rho, _ = spearmanr(actual, preds)
        print(f"{name}: Predicted Rank={' (Reversed)' if rho < 0 else ' (Correct)'}, rho={rho:.3f} {note}")

if __name__ == "__main__":
    simulate_rho_analysis()
