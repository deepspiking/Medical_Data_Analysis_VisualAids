import pandas as pd
import numpy as np
np.float = float
np.int = int
np.bool = bool
np.object = object
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr
import os

def main():
    os.makedirs('outputs/predictive_modeling', exist_ok=True)
    os.makedirs('outputs/pipeline_state', exist_ok=True)

    print("Defining X and y...")
    pheno_df = pd.read_csv('BRCA/BRCA_phenotype.txt', sep='\t', index_col='idx')
    rna_df = pd.read_csv('BRCA/BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt', sep='\t', index_col=0).T
    
    pheno_df.index = pheno_df.index.astype(str)
    rna_df.index = rna_df.index.astype(str)
    
    # y = High/Low ImmuneScore
    median_immune = pheno_df['ESTIMATE_ImmuneScore'].median()
    pheno_df['Target_y'] = (pheno_df['ESTIMATE_ImmuneScore'] > median_immune).astype(int)
    
    df_merged = rna_df.join(pheno_df['Target_y'], how='inner').dropna()
    X = df_merged.drop(columns=['Target_y'])
    y = df_merged['Target_y'].values
    
    # 1. Feature Selection (Correlation filter + LASSO)
    print("Filtering features by correlation with y...")
    corrs = []
    for col in X.columns:
        r, p = pearsonr(X[col], y)
        if p < 0.05:
            corrs.append((col, abs(r)))
    
    corrs.sort(key=lambda x: x[1], reverse=True)
    top_candidates = [x[0] for x in corrs[:300]]
    X_filtered = X[top_candidates]
    
    print("Standardizing and running LASSO...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_filtered)
    
    lasso_cv = LogisticRegressionCV(cv=5, penalty='l1', solver='liblinear', random_state=42, max_iter=1000)
    lasso_cv.fit(X_scaled, y)
    
    coefs = lasso_cv.coef_[0]
    selected_genes = X_filtered.columns[coefs != 0].tolist()
    
    if len(selected_genes) < 3 or len(selected_genes) > 20:
        top_indices = np.argsort(np.abs(coefs))[-10:]
        selected_genes = X_filtered.columns[top_indices].tolist()
        valid_coefs = coefs[top_indices]
    else:
        valid_coefs = coefs[coefs != 0]

    print(f"Selected {len(selected_genes)} genes: {selected_genes}")
    
    with open('outputs/pipeline_state/selected_features.txt', 'w') as f:
        f.write('\n'.join(selected_genes))
        
    df_merged[['Target_y'] + selected_genes].to_csv('outputs/pipeline_state/master_data.csv')
    
    plt.barh(selected_genes, valid_coefs)
    plt.title('LASSO Selected Biomarkers (Target: Immune Subtype)')
    plt.xlabel('Coefficient')
    plt.tight_layout()
    plt.savefig('outputs/predictive_modeling/lasso_path_plot.png')
    plt.close()

if __name__ == "__main__":
    main()
