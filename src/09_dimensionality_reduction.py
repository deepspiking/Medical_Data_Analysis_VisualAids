import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

def main():
    # Define file paths
    pheno_path = "BRCA/BRCA_phenotype.txt"
    rnaseq_path = "BRCA/BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"
    proteomics_path = "BRCA/BRCA_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
    
    out_dir = "outputs/advanced_plots"
    out_file = os.path.join(out_dir, "pca_tsne_plot.png")
    
    # Must do: Create outputs directory
    os.makedirs(out_dir, exist_ok=True)
    
    # Check which omics file is available
    omics_path = None
    if os.path.exists(rnaseq_path):
        omics_path = rnaseq_path
    elif os.path.exists(proteomics_path):
        omics_path = proteomics_path
        
    if not os.path.exists(pheno_path) or not omics_path:
        print(f"Warning: Ensure {pheno_path} and one of {rnaseq_path}/{proteomics_path} exist.")
        # Proceeding with mock data purely for demonstration if files are absent,
        # but in a real setting, we'd exit. We'll exit to conform to standard scripts.
        return
        
    try:
        # Load data (assuming tab-separated based on typical biological data)
        pheno = pd.read_csv(pheno_path, sep='\t', index_col=0)
        omics = pd.read_csv(omics_path, sep='\t', index_col=0)
    except Exception as e:
        print(f"Error loading data: {e}")
        return
        
    # Standardize orientation: Assume rows are features and columns are samples if features > samples
    if omics.shape[0] > omics.shape[1]:
        omics = omics.T
        
    # Align samples between phenotype and omics data
    common_samples = omics.index.intersection(pheno.index)
    if len(common_samples) == 0:
        print("Error: No overlapping samples between phenotype and omics data.")
        return
        
    omics = omics.loc[common_samples]
    pheno = pheno.loc[common_samples]
    
    # Must do: Handle missing values
    imputer = SimpleImputer(strategy='mean')
    omics_imputed = imputer.fit_transform(omics)
    
    # Scale data (PCA is sensitive to scaling)
    scaler = StandardScaler()
    omics_scaled = scaler.fit_transform(omics_imputed)
    
    # Dimensionality Reduction: PCA
    pca = PCA(n_components=2, random_state=42)
    pca_result = pca.fit_transform(omics_scaled)
    
    # Dimensionality Reduction: t-SNE
    # dynamically adjust perplexity for small sample sizes
    perplexity = min(30, max(5, len(common_samples) // 3))
    tsne = TSNE(n_components=2, perplexity=perplexity, random_state=42)
    tsne_result = tsne.fit_transform(omics_scaled)
    
    # Combine results for plotting
    plot_df = pd.DataFrame({
        'PCA1': pca_result[:, 0],
        'PCA2': pca_result[:, 1],
        'tSNE1': tsne_result[:, 0],
        'tSNE2': tsne_result[:, 1]
    }, index=common_samples)
    
    if 'ESTIMATE_ImmuneScore' in pheno.columns:
        median_immune = pheno['ESTIMATE_ImmuneScore'].median()
        plot_df['Immune_Subtype'] = np.where(pheno['ESTIMATE_ImmuneScore'] > median_immune, 'High Immune', 'Low Immune')
        target_col = 'Immune_Subtype'
    elif 'MSI_H' in pheno.columns:
        plot_df['MSI_Status'] = pheno['MSI_H'].fillna('Unknown')
        target_col = 'MSI_Status'
    else:
        target_col = None
        
    # Create a 1x2 subplot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # PCA Plot
    if target_col:
        sns.scatterplot(x='PCA1', y='PCA2', hue=target_col, data=plot_df, ax=axes[0], palette='Set1', alpha=0.8, s=60)
    else:
        sns.scatterplot(x='PCA1', y='PCA2', data=plot_df, ax=axes[0], alpha=0.8, s=60)
        
    axes[0].set_title('PCA of Omics Data', fontsize=14)
    axes[0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    axes[0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    
    # t-SNE Plot
    if target_col:
        sns.scatterplot(x='tSNE1', y='tSNE2', hue=target_col, data=plot_df, ax=axes[1], palette='Set1', alpha=0.8, s=60)
    else:
        sns.scatterplot(x='tSNE1', y='tSNE2', data=plot_df, ax=axes[1], alpha=0.8, s=60)
        
    axes[1].set_title('t-SNE of Omics Data', fontsize=14)
    
    # Save the output
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    print(f"Plot successfully saved to {out_file}")

if __name__ == "__main__":
    main()
