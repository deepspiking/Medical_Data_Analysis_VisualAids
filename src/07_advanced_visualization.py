import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def my_calibration_curve(y_true, y_prob, n_bins=10):
    bins = np.linspace(0., 1. + 1e-8, n_bins + 1)
    binids = np.digitize(y_prob, bins) - 1
    
    bin_sums = np.bincount(binids, weights=y_prob, minlength=len(bins))
    bin_true = np.bincount(binids, weights=y_true, minlength=len(bins))
    bin_total = np.bincount(binids, minlength=len(bins))
    
    nonzero = bin_total != 0
    prob_true = bin_true[nonzero] / bin_total[nonzero]
    prob_pred = bin_sums[nonzero] / bin_total[nonzero]
    
    return prob_true, prob_pred

def plot_nomogram_mock(output_path):
    """Generates a mock Nomogram template."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    variables = ['Age', 'Tumor Size', 'Grade', 'Total Points', 'Survival Probability']
    y_pos = np.arange(len(variables))
    
    for i, var in enumerate(variables):
        ax.plot([0, 100], [i, i], color='black')
        ax.text(-5, i, var, va='center', ha='right', fontsize=12)
        
        # Add some mock ticks
        for tick in range(0, 101, 10):
            ax.plot([tick, tick], [i - 0.1, i + 0.1], color='black')
            if tick % 20 == 0:
                ax.text(tick, i + 0.15, str(tick), va='bottom', ha='center', fontsize=8)

    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlim(-20, 110)
    ax.set_ylim(-1, len(variables))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    plt.title('Mock Nomogram Template', pad=20)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_calibration_curve_custom(y_true, y_prob, output_path):
    """Generates a Calibration Curve."""
    prob_true, prob_pred = my_calibration_curve(y_true, y_prob, n_bins=10)
    
    plt.figure(figsize=(8, 6))
    plt.plot(prob_pred, prob_true, marker='o', label='Model')
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray', label='Perfectly Calibrated')
    
    plt.xlabel('Mean Predicted Probability')
    plt.ylabel('Fraction of Positives')
    plt.title('Calibration Curve')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_dca(y_true, y_prob, output_path):
    """Generates a Decision Curve Analysis (DCA) Plot."""
    thresholds = np.linspace(0.01, 0.99, 100)
    net_benefit = []
    net_benefit_all = []
    
    N = len(y_true)
    
    for pt in thresholds:
        # Model net benefit
        tp = np.sum((y_prob >= pt) & (y_true == 1))
        fp = np.sum((y_prob >= pt) & (y_true == 0))
        nb = (tp / N) - (fp / N) * (pt / (1 - pt))
        net_benefit.append(nb)
        
        # Treat all net benefit
        tp_all = np.sum(y_true == 1)
        fp_all = np.sum(y_true == 0)
        nb_all = (tp_all / N) - (fp_all / N) * (pt / (1 - pt))
        net_benefit_all.append(nb_all)

    plt.figure(figsize=(8, 6))
    plt.plot(thresholds, net_benefit, label='Model', color='blue', linewidth=2)
    plt.plot(thresholds, net_benefit_all, label='Treat All', color='gray', linestyle='--')
    plt.plot(thresholds, np.zeros_like(thresholds), label='Treat None', color='black', linestyle=':')
    
    plt.ylim(-0.05, max(max(net_benefit), max(net_benefit_all)) + 0.05)
    plt.xlabel('Threshold Probability')
    plt.ylabel('Net Benefit')
    plt.title('Decision Curve Analysis (DCA)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_volcano(df, output_path):
    """Generates a Volcano Plot."""
    plt.figure(figsize=(8, 6))
    
    # Calculate -log10(p-value)
    df['neg_log_pval'] = -np.log10(df['p_value'])
    
    # Define colors based on thresholds
    p_thresh = -np.log10(0.05)
    fc_thresh = 1.0
    
    conditions = [
        (df['neg_log_pval'] > p_thresh) & (df['log2FC'] > fc_thresh),
        (df['neg_log_pval'] > p_thresh) & (df['log2FC'] < -fc_thresh)
    ]
    choices = ['Up-regulated', 'Down-regulated']
    df['Status'] = np.select(conditions, choices, default='Not Significant')
    
    colors = {'Up-regulated': 'red', 'Down-regulated': 'blue', 'Not Significant': 'gray'}
    
    sns.scatterplot(data=df, x='log2FC', y='neg_log_pval', hue='Status', palette=colors, alpha=0.7)
    
    plt.axvline(x=fc_thresh, color='k', linestyle='--', alpha=0.5)
    plt.axvline(x=-fc_thresh, color='k', linestyle='--', alpha=0.5)
    plt.axhline(y=p_thresh, color='k', linestyle='--', alpha=0.5)
    
    plt.xlabel('log2(Fold Change)')
    plt.ylabel('-log10(p-value)')
    plt.title('Volcano Plot')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_bubble(df, output_path):
    """Generates a Bubble Plot."""
    plt.figure(figsize=(10, 6))
    
    sns.scatterplot(data=df, x='Pathway_Enrichment', y='Significance', 
                    size='Gene_Count', sizes=(20, 500), hue='Category', 
                    alpha=0.7, palette='viridis')
    
    plt.xlabel('Pathway Enrichment Score')
    plt.ylabel('Significance (-log10 P)')
    plt.title('Bubble Plot of Pathway Enrichment')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_gradcam_mock(output_path):
    """Generates a Grad-CAM mock template."""
    plt.figure(figsize=(8, 8))
    
    # Create mock image (e.g., medical scan structure)
    base_image = np.zeros((100, 100))
    cv_x, cv_y = np.meshgrid(np.arange(100), np.arange(100))
    base_image = np.exp(-((cv_x - 50)**2 + (cv_y - 50)**2) / (2 * 20**2))
    
    # Create mock heatmap
    heatmap = np.exp(-((cv_x - 60)**2 + (cv_y - 40)**2) / (2 * 10**2))
    
    # Plot base image in grayscale
    plt.imshow(base_image, cmap='gray')
    
    # Overlay heatmap with jet colormap and transparency
    plt.imshow(heatmap, cmap='jet', alpha=0.5)
    
    plt.colorbar(label='Activation Intensity')
    plt.title('Grad-CAM Mock Template (Activation Map)')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

if __name__ == '__main__':
    out_dir = 'outputs/advanced_plots'
    os.makedirs(out_dir, exist_ok=True)
    
    print("Generating Nomogram mock...")
    plot_nomogram_mock(f'{out_dir}/nomogram_mock.png')
    
    print("Generating Calibration Curve...")
    # Synthetic data for calibration and DCA
    np.random.seed(42)
    y_true = np.random.binomial(1, 0.3, 1000)
    y_prob = y_true * np.random.beta(5, 2, 1000) + (1 - y_true) * np.random.beta(2, 5, 1000)
    y_prob = np.clip(y_prob, 0, 1)
    plot_calibration_curve_custom(y_true, y_prob, f'{out_dir}/calibration_curve.png')
    
    print("Generating DCA Plot...")
    plot_dca(y_true, y_prob, f'{out_dir}/dca_plot.png')
    
    print("Generating Volcano Plot...")
    # Synthetic data for Volcano plot
    genes_df = pd.DataFrame({
        'Gene': [f'Gene_{i}' for i in range(500)],
        'log2FC': np.random.normal(0, 1.5, 500),
        'p_value': np.random.uniform(0.0001, 1.0, 500)
    })
    # Make some genes clearly significant
    genes_df.loc[0:20, 'log2FC'] = np.random.normal(2.5, 0.5, 21)
    genes_df.loc[0:20, 'p_value'] = np.random.uniform(0.00001, 0.01, 21)
    genes_df.loc[21:40, 'log2FC'] = np.random.normal(-2.5, 0.5, 20)
    genes_df.loc[21:40, 'p_value'] = np.random.uniform(0.00001, 0.01, 20)
    plot_volcano(genes_df, f'{out_dir}/volcano_plot.png')
    
    print("Generating Bubble Plot...")
    # Synthetic data for Bubble plot
    bubble_df = pd.DataFrame({
        'Pathway_Enrichment': np.random.uniform(0.5, 5.0, 30),
        'Significance': np.random.uniform(1.0, 10.0, 30),
        'Gene_Count': np.random.randint(10, 200, 30),
        'Category': np.random.choice(['Immune', 'Metabolic', 'Signaling', 'Cell Cycle'], 30)
    })
    plot_bubble(bubble_df, f'{out_dir}/bubble_plot.png')
    
    print("Generating Grad-CAM Mock...")
    plot_gradcam_mock(f'{out_dir}/gradcam_mock.png')
    
    print(f"All plots generated successfully in '{out_dir}' directory.")
