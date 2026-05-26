import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    # Create output directory
    out_dir = "outputs/advanced_plots"
    os.makedirs(out_dir, exist_ok=True)

    # Mock data for pathways
    data = {
        'Pathway': [
            'p53 signaling pathway',
            'Apoptosis',
            'PI3K-Akt signaling pathway',
            'mTOR signaling pathway',
            'Cell cycle',
            'Notch signaling pathway',
            'MAPK signaling pathway',
            'Wnt signaling pathway'
        ],
        'Fold_Enrichment': [4.5, 4.2, 3.8, 3.5, 2.9, 2.4, 2.1, 1.8],
        'p_value': [5e-6, 1e-5, 2e-4, 1e-4, 1e-3, 0.005, 5e-3, 0.01]
    }
    
    df = pd.DataFrame(data)
    # Calculate -log10(p_value)
    df['-log10(p_value)'] = -np.log10(df['p_value'])
    
    # Sort by -log10(p_value) for better visualization (highest on top)
    df = df.sort_values('-log10(p_value)', ascending=False)
    
    # Set the style
    sns.set_theme(style="whitegrid")
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Use a color palette that maps to Fold Enrichment
    norm = plt.Normalize(df['Fold_Enrichment'].min(), df['Fold_Enrichment'].max())
    cmap = plt.get_cmap("YlOrRd")
    colors = [cmap(norm(val)) for val in df['Fold_Enrichment']]
    
    sns.barplot(
        x='-log10(p_value)', 
        y='Pathway', 
        data=df,
        palette=colors,
        edgecolor='black',
        ax=ax
    )
    
    # Add colorbar for Fold Enrichment
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Fold Enrichment', rotation=270, labelpad=15)
    
    # Set labels and title
    ax.set_xlabel('-log10(p-value)', fontsize=12)
    ax.set_ylabel('')
    ax.set_title('Pathway Enrichment Analysis', fontsize=14, pad=15)
    
    plt.tight_layout()
    
    # Save the plot
    out_path = os.path.join(out_dir, "gsea_barplot.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"Saved mock GSEA barplot to {out_path}")

if __name__ == "__main__":
    main()
