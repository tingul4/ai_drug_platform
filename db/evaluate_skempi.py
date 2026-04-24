
import sqlite3
import numpy as np
from scipy.stats import spearmanr, pearsonr
from engine.analyzer import SKEMPIModel
import matplotlib.pyplot as plt

def evaluate():
    conn = sqlite3.connect("db/skempi.db")
    conn.row_factory = sqlite3.Row
    model = SKEMPIModel(conn)
    
    # Fetch single-point mutations with experimental ddG
    cur = conn.execute("""
        SELECT mutation_clean, ddG_kcal, location
        FROM skempi_mutations
        WHERE ddG_kcal IS NOT NULL
          AND mutation_clean NOT LIKE '%,%' -- Single point only
          AND length(mutation_clean) >= 3
    """)
    
    y_true = []
    y_pred = []
    
    for row in cur:
        mut = row['mutation_clean'].strip()
        wt = mut[0]
        mut_aa = mut[-1]
        
        # Simple prediction based on the model logic
        # SKEMPIModel.predict_ddG uses these parameters
        ddG_pred, confidence, n = model.predict_ddG(wt, mut_aa, is_interface=True)
        
        y_true.append(row['ddG_kcal'])
        y_pred.append(ddG_pred)
    
    conn.close()
    
    if not y_true:
        print("No data found for evaluation.")
        return

    rho, p_val = spearmanr(y_true, y_pred)
    r, _ = pearsonr(y_true, y_pred)
    
    print(f"\n=== SKEMPI Global Performance (n={len(y_true)}) ===")
    print(f"Spearman rho: {rho:.3f} (p={p_val:.2e})")
    print(f"Pearson r:    {r:.3f}")
    print(f"Mean Abs Error: {np.mean(np.abs(np.array(y_true) - np.array(y_pred))):.3f} kcal/mol")
    
    # Optional: logic check for PAI-1 points
    print("\n=== Mechanism Check (Directional Accuracy) ===")
    correct_direction = sum(1 for t, p in zip(y_true, y_pred) if (t > 0.5 and p > 0) or (t < -0.5 and p < 0))
    significant_muts = sum(1 for t in y_true if abs(t) > 0.5)
    print(f"Correct Direction on significant mutations (|ddG|>0.5): {correct_direction}/{significant_muts} ({correct_direction/significant_muts*100:.1f}%)")

if __name__ == "__main__":
    evaluate()
