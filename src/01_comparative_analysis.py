import os
import pandas as pd
import numpy as np
from scipy import stats
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def load_data(phenotype_path, rnaseq_path):
    """
    Load phenotype and RNAseq data, and merge them by sample ID.
    Handles potential file formats correctly and transposes RNAseq data if needed.
    """
    if not os.path.exists(phenotype_path):
        print(f"Warning: {phenotype_path} not found.")
        phenotype_df = pd.DataFrame()
    else:
        # Load phenotype data
        phenotype_df = pd.read_csv(phenotype_path, sep='\t', index_col=0)

    if not os.path.exists(rnaseq_path):
        print(f"Warning: {rnaseq_path} not found.")
        rnaseq_df = pd.DataFrame()
    else:
        # Load RNAseq data
        rnaseq_df = pd.read_csv(rnaseq_path, sep='\t', index_col=0)
        # Transpose to have samples as rows and genes as columns
        rnaseq_df = rnaseq_df.T

    # Merge on index (sample ID)
    if not phenotype_df.empty and not rnaseq_df.empty:
        # Keep inner join of samples present in both
        merged_df = phenotype_df.join(rnaseq_df, how='inner')
    else:
        merged_df = pd.DataFrame()
        
    return merged_df, phenotype_df, rnaseq_df

def run_t_test(group1, group2):
    """
    Perform Student's t-test for independent samples.
    """
    # Drop NaNs to handle missing data gracefully
    group1 = group1.dropna()
    group2 = group2.dropna()
    if len(group1) == 0 or len(group2) == 0:
        return np.nan, np.nan
    # Welch's t-test is generally safer as it does not assume equal variances
    t_stat, p_val = stats.ttest_ind(group1, group2, equal_var=False)
    return t_stat, p_val

def run_mann_whitney(group1, group2):
    """
    Perform Mann-Whitney U test (non-parametric equivalent of independent t-test).
    """
    group1 = group1.dropna()
    group2 = group2.dropna()
    if len(group1) == 0 or len(group2) == 0:
        return np.nan, np.nan
    u_stat, p_val = stats.mannwhitneyu(group1, group2, alternative='two-sided')
    return u_stat, p_val

def run_anova(*groups):
    """
    Perform One-Way ANOVA.
    """
    clean_groups = [g.dropna() for g in groups]
    clean_groups = [g for g in clean_groups if len(g) > 0]
    if len(clean_groups) < 2:
        return np.nan, np.nan
    f_stat, p_val = stats.f_oneway(*clean_groups)
    return f_stat, p_val

def run_kruskal_wallis(*groups):
    """
    Perform Kruskal-Wallis H test (non-parametric equivalent of One-Way ANOVA).
    """
    clean_groups = [g.dropna() for g in groups]
    clean_groups = [g for g in clean_groups if len(g) > 0]
    if len(clean_groups) < 2:
        return np.nan, np.nan
    h_stat, p_val = stats.kruskal(*clean_groups)
    return h_stat, p_val

def main():
    parser = argparse.ArgumentParser(description="Run comparative analysis on BRCA data")
    parser.add_argument('--pheno', default='BRCA/BRCA_phenotype.txt', help='Path to phenotype data')
    parser.add_argument('--rna', default='BRCA/BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt', help='Path to RNAseq data')
    
    args = parser.parse_args()
    
    merged_df, pheno, rna = load_data(args.pheno, args.rna)
    
    if merged_df.empty:
        print("Merged data is empty or missing. Please ensure the files exist.")
        return
        
    print(f"Data loaded successfully. Total samples: {len(merged_df)}")
    
    # Example usage on the first gene if available
    if len(rna.columns) > 0 and len(pheno.columns) > 0:
        gene = rna.columns[0]
        # Find a suitable categorical column with a few unique values (e.g., stage, subtype)
        cat_col = None
        for col in pheno.columns:
            unique_vals = pheno[col].nunique()
            if 1 < unique_vals <= 5:
                cat_col = col
                break
                
        if cat_col:
            print(f"\nPerforming example analysis on gene '{gene}' grouped by phenotype '{cat_col}'")
            
            plot_df = merged_df[[cat_col, gene]].dropna()
            
            unique_groups = plot_df[cat_col].unique()
            groups = [plot_df[plot_df[cat_col] == g][gene] for g in unique_groups]
            
            if len(groups) == 2:
                t_stat, p_val = run_t_test(groups[0], groups[1])
                print(f"T-test: t={t_stat:.4f}, p={p_val:.4g}")
                
                u_stat, u_p_val = run_mann_whitney(groups[0], groups[1])
                print(f"Mann-Whitney U: U={u_stat:.4f}, p={u_p_val:.4g}")

                out_dir = 'outputs/comparative_analysis'
                os.makedirs(out_dir, exist_ok=True)
                
                plt.figure(figsize=(8, 6))
                ax = sns.boxplot(x=cat_col, y=gene, data=plot_df, order=unique_groups)
                
                x1, x2 = 0, 1
                y = plot_df[gene].max() + np.abs(plot_df[gene].max() * 0.05) + 0.1
                h = np.abs(plot_df[gene].max() * 0.02) + 0.1
                col = 'k'
                
                plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
                plt.text((x1+x2)*.5, y+h, f"p = {p_val:.2e}", ha='center', va='bottom', color=col)
                
                plt.title(f"{gene} Expression by {cat_col}")
                plt.ylabel(f"{gene} Expression")
                plt.tight_layout()
                
                out_path = os.path.join(out_dir, 'boxplot_pvalue.png')
                plt.savefig(out_path)
                plt.close()
                print(f"Boxplot saved to {out_path}")
                
            elif len(groups) > 2:
                f_stat, p_val = run_anova(*groups)
                print(f"ANOVA: F={f_stat:.4f}, p={p_val:.4g}")
                
                h_stat, p_val = run_kruskal_wallis(*groups)
                print(f"Kruskal-Wallis: H={h_stat:.4f}, p={p_val:.4g}")
        else:
            print("No suitable categorical phenotype column found for example analysis.")
    else:
        print("Missing genes or phenotype columns for analysis.")

if __name__ == "__main__":
    main()
