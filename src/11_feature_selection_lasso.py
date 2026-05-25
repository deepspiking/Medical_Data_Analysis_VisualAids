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
    plt.close()
    print(f"Plot saved to {output_file}")
    
    print("Extracting LASSO-selected features...")
    selected_indices = np.where(lasso_cv.coef_ != 0)[0]
    selected_genes = top_500_genes[selected_indices]
    
    if len(selected_genes) < 2:
        print("LASSO dropped all features at optimal alpha. Selecting top 10 features from a looser alpha for demonstration...")
        idx = len(lasso_cv.alphas_) // 2
        path_coefs = coefs_path[:, idx]
        top_indices = np.argsort(np.abs(path_coefs))[-10:]
        selected_genes = top_500_genes[top_indices]
        
    print(f"Number of features selected for correlation: {len(selected_genes)}")
    
    if len(selected_genes) > 1:
        print("Generating correlation heatmap for selected features...")
        import seaborn as sns
        from scipy import stats
        
        selected_data = df_merged[list(selected_genes) + ['OS_days']]
        corr_matrix = selected_data.corr(method='spearman')
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr_matrix, annot=True, fmt=".2f", cmap='coolwarm', center=0, 
                    square=True, linewidths=.5, cbar_kws={"shrink": .8})
        plt.title('Correlation Heatmap of LASSO-Selected Features & OS_days', fontsize=14, pad=20)
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        
        corr_output_file = 'outputs/predictive_modeling/lasso_selected_features_corr.png'
        plt.savefig(corr_output_file, bbox_inches='tight', dpi=300)
        plt.close()
        print(f"Correlation heatmap of selected features saved to {corr_output_file}")
    else:
        print("Not enough features selected to generate a correlation heatmap.")

if __name__ == "__main__":
    main()
