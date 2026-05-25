import pandas as pd
import numpy as np
np.float = float
np.int = int
np.bool = bool
import matplotlib.pyplot as plt
from sklearn.linear_model import LassoCV, lasso_path
from sklearn.preprocessing import StandardScaler
import os

def main():
    # Ensure output directory exists
    os.makedirs('outputs/predictive_modeling', exist_ok=True)

    print("Loading RNAseq data...")
    # Read RNAseq data
    # format: genes as rows, samples as columns. First column is 'idx'
    df_rna = pd.read_csv('BRCA/BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt', sep='\t', index_col=0)
    
    print("Transposing and filtering top 500 variable genes...")
    # Transpose so samples are rows and genes are columns
    df_rna = df_rna.T
    
    # Calculate variance of each gene and take top 500
    variances = df_rna.var()
    top_500_genes = variances.nlargest(500).index
    df_rna_filtered = df_rna[top_500_genes]

    print("Loading survival data...")
    df_surv = pd.read_csv('BRCA/BRCA_survival.txt', sep='\t', index_col=0)

    # Merge data
    print("Merging data...")
    # df_surv has index as case_id
    # df_rna has index as case_id
    df_merged = df_rna_filtered.join(df_surv['OS_days'], how='inner')
    
    # Drop NaNs in target
    df_merged = df_merged.dropna(subset=['OS_days'])
    
    X = df_merged.drop(columns=['OS_days']).values
    y = df_merged['OS_days'].values
    
    print("Standardizing features...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    print("Computing LASSO path...")
    alphas_path, coefs_path, _ = lasso_path(X_scaled, y)
    
    print("Running LassoCV to find optimal alpha...")
    lasso_cv = LassoCV(cv=5, random_state=42)
    lasso_cv.fit(X_scaled, y)
    
    # Plotting
    print("Plotting LASSO path...")
    plt.figure(figsize=(10, 6))
    
    # Plot paths
    # coefs_path is (n_features, n_alphas)
    log_alphas = np.log10(alphas_path)
    for coef in coefs_path:
        plt.plot(log_alphas, coef)
        
    plt.axvline(np.log10(lasso_cv.alpha_), linestyle='--', color='k', label='Optimal Alpha (CV)')
    
    plt.xlabel('Log10(Alpha)')
    plt.ylabel('Coefficients')
    plt.title('LASSO Coefficient Path')
    plt.legend()
    plt.grid(True)
    
    output_file = 'outputs/predictive_modeling/lasso_path_plot.png'
    plt.savefig(output_file, bbox_inches='tight')
    print(f"Plot saved to {output_file}")

if __name__ == "__main__":
    main()
