import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

def main():
    # File paths
    pheno_file = "BRCA/BRCA_phenotype.txt"
    rna_file = "BRCA/BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"
    out_dir = "outputs/advanced_plots"
    out_file = os.path.join(out_dir, "hierarchical_clustering.png")

    os.makedirs(out_dir, exist_ok=True)

    print("Loading phenotype data...")
    pheno_df = pd.read_csv(pheno_file, sep='\t', index_col='idx')
    
    # Create a categorical group based on ESTIMATE_StromalScore
    median_score = pheno_df['ESTIMATE_StromalScore'].median()
    pheno_df['Stromal_Category'] = np.where(pheno_df['ESTIMATE_StromalScore'] >= median_score, 'High', 'Low')

    print("Loading RNAseq data...")
    rna_df = pd.read_csv(rna_file, sep='\t', index_col='idx')

    # Calculate variance for each gene across samples
    print("Selecting top 50 highly variable genes...")
    gene_vars = rna_df.var(axis=1)
    top_50_genes = gene_vars.nlargest(50).index

    # Filter RNAseq data to top 50 genes and transpose so samples are rows
    rna_top = rna_df.loc[top_50_genes].T

    # Align with phenotype data
    common_samples = rna_top.index.intersection(pheno_df.index)
    rna_top = rna_top.loc[common_samples]
    pheno_aligned = pheno_df.loc[common_samples]

    # Map colors to Stromal_Category
    categories = pheno_aligned['Stromal_Category'].unique()
    # Simple color palette mapping
    palette = sns.color_palette("Set2", len(categories))
    lut = dict(zip(categories, palette))
    row_colors = pheno_aligned['Stromal_Category'].map(lut)

    print("Generating clustermap...")
    # Generate hierarchical clustering heatmap
    g = sns.clustermap(
        rna_top,
        cmap="viridis",
        row_colors=row_colors,
        xticklabels=True,  # 50 genes are readable
        yticklabels=False,  # Hide sample labels
        figsize=(12, 10),
        z_score=1,  # Standardize across columns (genes) using z-score is often better for heatmaps
        method='ward'
    )

    # Add legend for the row colors
    for label in categories:
        g.ax_col_dendrogram.bar(0, 0, color=lut[label], label=label, linewidth=0)
    g.ax_col_dendrogram.legend(title="Stromal Category", loc="center", ncol=2, bbox_to_anchor=(0.5, 1.2))

    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    print(f"Saved plot to {out_file}")

if __name__ == "__main__":
    main()
