import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import export_graphviz
import requests

def main():
    out_dir = 'outputs/predictive_modeling'
    os.makedirs(out_dir, exist_ok=True)
    
    # Read pipeline state
    feat_file = 'outputs/pipeline_state/selected_features.txt'
    if not os.path.exists(feat_file):
        print(f"Run script 11 first. {feat_file} not found.")
        return
    with open(feat_file, 'r') as f:
        selected_genes = [line.strip() for line in f if line.strip()]

    print("Loading data...")
    rna_df = pd.read_csv('BRCA/BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt', sep='\t', index_col=0).T
    surv_df = pd.read_csv('BRCA/BRCA_survival.txt', sep='\t', index_col=0)
    
    merged_df = rna_df[selected_genes].join(surv_df, how='inner').dropna(subset=['OS_event'])
    
    X = merged_df[selected_genes]
    y = merged_df['OS_event'].astype(int)
    
    print("Training Decision Tree on LASSO features...")
    dt_model = DecisionTreeClassifier(max_depth=3, random_state=42, class_weight='balanced')
    dt_model.fit(X, y)
    
    # Save Risk Scores for Survival Analysis (Step 03)
    preds = dt_model.predict(X)
    risk_df = pd.DataFrame({'Risk_Group': np.where(preds == 1, 'High Risk', 'Low Risk')}, index=X.index)
    risk_df.to_csv('outputs/pipeline_state/patient_risk.csv')
    print("Saved patient risk groups to pipeline_state/patient_risk.csv")
    
    # Export Graphviz
    dot_data = export_graphviz(dt_model, out_file=None, 
                               feature_names=selected_genes, 
                               class_names=['Alive', 'Deceased'], 
                               filled=True, rounded=True, special_characters=True)
    try:
        resp = requests.post("https://quickchart.io/graphviz", json={"graph": dot_data}, timeout=10)
        if resp.status_code == 200:
            with open(os.path.join(out_dir, 'decision_tree_viz.svg'), 'wb') as f:
                f.write(resp.content)
            print("Saved Decision Tree SVG.")
    except Exception as e:
        print(f"Graphviz API failed: {e}")

if __name__ == "__main__":
    main()
