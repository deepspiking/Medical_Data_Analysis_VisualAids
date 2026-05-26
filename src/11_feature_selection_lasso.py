import pandas as pd
import numpy as np
np.float = float
np.int = int
np.bool = bool
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
import os

def main():
    os.makedirs('outputs/predictive_modeling', exist_ok=True)
    os.makedirs('outputs/pipeline_state', exist_ok=True)

    print("Loading RNAseq data...")
    df_rna = pd.read_csv('BRCA/BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt', sep='\t', index_col=0).T
    
    # Take top 500 variable genes to speed up LASSO
    top_genes = df_rna.var().nlargest(500).index
    df_rna_filtered = df_rna[top_genes]

    print("Loading survival data...")
    df_surv = pd.read_csv('BRCA/BRCA_survival.txt', sep='\t', index_col=0)
    
    # Merge
    df_merged = df_rna_filtered.join(df_surv['OS_event'], how='inner').dropna(subset=['OS_event'])
    
    X = df_merged.drop(columns=['OS_event'])
    y = df_merged['OS_event'].values
    
    print("Standardizing features & Running LASSO (Logistic L1)...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Use LogisticRegressionCV with L1 penalty for classification feature selection
    lasso_cv = LogisticRegressionCV(cv=5, penalty='l1', solver='liblinear', random_state=42, class_weight='balanced', max_iter=1000)
    lasso_cv.fit(X_scaled, y)
    
    # Extract selected features
    coefs = lasso_cv.coef_[0]
    selected_indices = np.where(coefs != 0)[0]
    selected_genes = X.columns[selected_indices].tolist()
    
    # Fallback if too few/many genes selected
    if len(selected_genes) < 3 or len(selected_genes) > 15:
        print(f"LASSO selected {len(selected_genes)} genes. Falling back to top 8 by absolute coefficient for stable pipeline...")
        top_indices = np.argsort(np.abs(coefs))[-8:]
        selected_genes = X.columns[top_indices].tolist()
        
    print(f"\nFinal Selected Features for downstream pipeline: {selected_genes}")
    
    # Save to state file for other scripts
    with open('outputs/pipeline_state/selected_features.txt', 'w') as f:
        f.write('\n'.join(selected_genes))
        
    # Dummy plot for compatibility
    plt.figure(figsize=(8, 5))
    plt.barh(selected_genes, coefs[np.isin(X.columns, selected_genes)])
    plt.title('LASSO Selected Features & Coefficients')
    plt.xlabel('Coefficient Value')
    plt.tight_layout()
    plt.savefig('outputs/predictive_modeling/lasso_path_plot.png')
    plt.close()
    
if __name__ == "__main__":
    main()
