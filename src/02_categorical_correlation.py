import os
import pandas as pd
import numpy as np
np.float = float
np.int = int
np.bool = bool
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def main():
    out_dir = 'outputs/categorical_correlation'
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
    
    rna_df.index = rna_df.index.astype(str)
    surv_df.index = surv_df.index.astype(str)
    
    merged_df = rna_df.join(surv_df, how='inner')
    
    print("\n--- Correlation Analysis of LASSO Features ---")
    corr_cols = ['OS_event', 'OS_days'] + selected_genes
    corr_data = merged_df[corr_cols].apply(pd.to_numeric, errors='coerce')
    
    corr_matrix = corr_data.corr(method='spearman')
    pval_matrix = pd.DataFrame(np.ones_like(corr_matrix.values), columns=corr_cols, index=corr_cols)
    for i in corr_cols:
        for j in corr_cols:
            if i != j:
                mask = corr_data[i].notna() & corr_data[j].notna()
                if mask.sum() > 2:
                    _, p_val = stats.spearmanr(corr_data.loc[mask, i], corr_data.loc[mask, j])
                    pval_matrix.loc[i, j] = p_val
                
    annot_matrix = pd.DataFrame('', index=corr_matrix.index, columns=corr_matrix.columns)
    for i in corr_cols:
        for j in corr_cols:
            r = corr_matrix.loc[i, j]
            p = pval_matrix.loc[i, j]
            stars = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            annot_matrix.loc[i, j] = '1.0' if i == j else f"{r:.2f}\n{stars}"
                
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr_matrix, annot=annot_matrix.values, fmt='', cmap='coolwarm', center=0, 
                square=True, linewidths=.5, cbar_kws={"shrink": .8})
    plt.title("Correlation Matrix of LASSO-Selected Features & Survival", fontsize=16, pad=20)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'correlation_heatmap.png'), dpi=300)
    plt.close()
    print("Saved correlation heatmap.")
    
    # Save a dummy categorical plot just to not break the markdown image links
    plt.figure(figsize=(6,4))
    plt.text(0.5, 0.5, "Replaced by unified pipeline", ha='center')
    plt.savefig(os.path.join(out_dir, 'categorical_barplot.png'))
    plt.close()

if __name__ == '__main__':
    main()
