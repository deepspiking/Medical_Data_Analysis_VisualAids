import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def main():
    # 1. Load data
    print("Loading phenotype data...")
    pheno_df = pd.read_csv('BRCA/BRCA_phenotype.txt', sep='\t', index_col='idx')
    
    print("Loading proteomics data...")
    prot_df = pd.read_csv('BRCA/BRCA_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt', sep='\t', index_col='idx')
    
    # Transpose proteomics data so samples are rows and genes are columns
    prot_df = prot_df.T
    
    # 2. Merge data on sample index
    print("Merging data...")
    # Ensure index types match (sometimes they are parsed differently)
    pheno_df.index = pheno_df.index.astype(str)
    prot_df.index = prot_df.index.astype(str)
    merged_df = pheno_df.join(prot_df, how='inner')
    print(f"Merged dataframe shape: {merged_df.shape}")
    
    # 3. Categorical Analysis
    print("\n--- Categorical Analysis ---")
    
    # Create categorical features from continuous ones to guarantee variance for tests
    print("Creating binary categories from continuous scores (High > Median, Low <= Median)...")
    merged_df['Stromal_Category'] = np.where(merged_df['ESTIMATE_StromalScore'] > merged_df['ESTIMATE_StromalScore'].median(), 'High', 'Low')
    merged_df['Immune_Category'] = np.where(merged_df['ESTIMATE_ImmuneScore'] > merged_df['ESTIMATE_ImmuneScore'].median(), 'High', 'Low')
    
    contingency_table = pd.crosstab(merged_df['Stromal_Category'], merged_df['Immune_Category'])
    print("\nCross-tabulation (Stromal vs Immune Categories):")
    print(contingency_table)
    
    # Chi-Square Test
    try:
        chi2, p_chi2, dof, expected = stats.chi2_contingency(contingency_table)
        print(f"\nChi-Square Test: chi2={chi2:.4f}, p-value={p_chi2:.4e}")
    except Exception as e:
        print(f"Chi-Square Test failed: {e}")
        
    # Fisher's Exact Test
    if contingency_table.shape == (2, 2):
        try:
            res = stats.fisher_exact(contingency_table)
            print(f"Fisher's Exact Test: odds_ratio={res[0]:.4f}, p-value={res[1]:.4e}")
        except Exception as e:
            print(f"Fisher's Exact Test failed: {e}")
    else:
        print("Fisher's Exact Test requires a 2x2 table.")

    # 4. Correlation Analysis
    print("\n--- Correlation Analysis ---")
    # Let's correlate the first gene with a phenotype score
    gene = prot_df.columns[2] # Picking the third gene just in case the first has many NAs
    pheno_feat = 'ESTIMATE_StromalScore'
    
    if gene in merged_df.columns and pheno_feat in merged_df.columns:
        valid_corr = merged_df[[gene, pheno_feat]].dropna()
        valid_corr[gene] = pd.to_numeric(valid_corr[gene], errors='coerce')
        valid_corr = valid_corr.dropna()
        
        print(f"Calculating correlations between {gene} and {pheno_feat} (N={len(valid_corr)})...")
        if len(valid_corr) > 2:
            x = valid_corr[gene]
            y = valid_corr[pheno_feat]
            
            # Pearson
            r_pearson, p_pearson = stats.pearsonr(x, y)
            print(f"Pearson Correlation: r={r_pearson:.4f}, p-value={p_pearson:.4e}")
            
            # Spearman
            r_spearman, p_spearman = stats.spearmanr(x, y)
            print(f"Spearman Correlation: rho={r_spearman:.4f}, p-value={p_spearman:.4e}")
            
            # Plotting
            out_dir = 'outputs/categorical_correlation'
            os.makedirs(out_dir, exist_ok=True)
            plt.figure(figsize=(8, 6))
            sns.regplot(x=x, y=y)
            plt.title(f"Correlation between {gene} and {pheno_feat}")
            plt.xlabel(gene)
            plt.ylabel(pheno_feat)
            plt.text(0.05, 0.95, f"Pearson r = {r_pearson:.2f}\np = {p_pearson:.2e}", 
                     transform=plt.gca().transAxes, fontsize=12, verticalalignment='top',
                     bbox=dict(boxstyle="round,pad=0.3", alpha=0.5, facecolor="white"))
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, 'scatter_corr.png'))
            plt.close()
            print(f"Saved scatter plot to {os.path.join(out_dir, 'scatter_corr.png')}")
        else:
            print("Not enough valid data points for correlation.")
    else:
        print("Specified columns for correlation not found in merged data.")

if __name__ == '__main__':
    main()
