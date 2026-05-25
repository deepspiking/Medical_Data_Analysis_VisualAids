import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def main():
    print("Loading phenotype data...")
    pheno_df = pd.read_csv('BRCA/BRCA_phenotype.txt', sep='\t', index_col='idx')
    
    print("Loading proteomics data...")
    prot_df = pd.read_csv('BRCA/BRCA_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt', sep='\t', index_col='idx')
    prot_df = prot_df.T
    
    pheno_df.index = pheno_df.index.astype(str)
    prot_df.index = prot_df.index.astype(str)
    merged_df = pheno_df.join(prot_df, how='inner')
    
    out_dir = 'outputs/categorical_correlation'
    os.makedirs(out_dir, exist_ok=True)
    
    # ---------------------------------------------------------
    # 1. Categorical Analysis Visualization (Stacked Bar Chart)
    # ---------------------------------------------------------
    print("\n--- Categorical Analysis ---")
    merged_df['Stromal_Category'] = np.where(merged_df['ESTIMATE_StromalScore'] > merged_df['ESTIMATE_StromalScore'].median(), 'High Stromal', 'Low Stromal')
    merged_df['Immune_Category'] = np.where(merged_df['ESTIMATE_ImmuneScore'] > merged_df['ESTIMATE_ImmuneScore'].median(), 'High Immune', 'Low Immune')
    
    contingency_table = pd.crosstab(merged_df['Stromal_Category'], merged_df['Immune_Category'])
    
    chi2, p_chi2, dof, expected = stats.chi2_contingency(contingency_table)
    
    # Convert to percentages for 100% stacked bar plot
    crosstab_pct = contingency_table.div(contingency_table.sum(1), axis=0) * 100
    
    fig, ax = plt.subplots(figsize=(8, 6))
    crosstab_pct.plot(kind='bar', stacked=True, color=['#1f77b4', '#ff7f0e'], ax=ax)
    
    plt.title("Categorical Distribution: Stromal vs Immune Category", fontsize=14, pad=20)
    plt.ylabel("Percentage (%)", fontsize=12)
    plt.xlabel("Stromal Category", fontsize=12)
    plt.xticks(rotation=0)
    plt.legend(title="Immune Category", bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Add p-value text
    p_text = f"Chi-square p-value: {p_chi2:.4f}"
    if p_chi2 < 0.05:
        p_text += " (*)"
    plt.text(0.5, 1.05, p_text, ha='center', va='center', transform=ax.transAxes, fontsize=12, fontweight='bold', color='darkred')
    
    plt.tight_layout()
    cat_out = os.path.join(out_dir, 'categorical_barplot.png')
    plt.savefig(cat_out, dpi=300)
    plt.close()
    print(f"Saved categorical bar plot to {cat_out}")

    # ---------------------------------------------------------
    # 2. Correlation Matrix Heatmap
    # ---------------------------------------------------------
    print("\n--- Correlation Analysis ---")
    # Select 8 highly variable genes and 2 phenotype scores
    top_genes = prot_df.var().sort_values(ascending=False).head(8).index.tolist()
    pheno_feats = ['ESTIMATE_StromalScore', 'ESTIMATE_ImmuneScore']
    corr_cols = pheno_feats + top_genes
    
    corr_data = merged_df[corr_cols].apply(pd.to_numeric, errors='coerce').dropna()
    
    # Calculate correlation matrix (Spearman for robustness)
    corr_matrix = corr_data.corr(method='spearman')
    
    # Calculate p-values to add asterisks
    pval_matrix = pd.DataFrame(np.ones_like(corr_matrix.values), columns=corr_cols, index=corr_cols)
    for i in corr_cols:
        for j in corr_cols:
            if i != j:
                _, p_val = stats.spearmanr(corr_data[i], corr_data[j])
                pval_matrix.loc[i, j] = p_val
                
    # Format annotations: r value + * if significant
    annot_matrix = pd.DataFrame('', index=corr_matrix.index, columns=corr_matrix.columns)
    for i in corr_cols:
        for j in corr_cols:
            r = corr_matrix.loc[i, j]
            p = pval_matrix.loc[i, j]
            stars = ''
            if p < 0.001: stars = '***'
            elif p < 0.01: stars = '**'
            elif p < 0.05: stars = '*'
            
            if i == j:
                annot_matrix.loc[i, j] = '1.0'
            else:
                annot_matrix.loc[i, j] = f"{r:.2f}\n{stars}"
                
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr_matrix, annot=annot_matrix.values, fmt='', cmap='coolwarm', vmin=-1, vmax=1, 
                square=True, linewidths=.5, cbar_kws={"shrink": .8})
    plt.title("Correlation Matrix Heatmap (Spearman rho)", fontsize=16, pad=20)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    
    plt.tight_layout()
    corr_out = os.path.join(out_dir, 'correlation_heatmap.png')
    plt.savefig(corr_out, dpi=300)
    plt.close()
    print(f"Saved correlation heatmap to {corr_out}")

if __name__ == '__main__':
    main()
