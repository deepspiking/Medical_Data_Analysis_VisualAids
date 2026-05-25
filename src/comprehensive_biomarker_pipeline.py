import pandas as pd
import numpy as np
np.float = float
np.int = int
np.bool = bool
np.object = object

import matplotlib.pyplot as plt
import seaborn as sns
import os
import requests

from scipy.stats import pearsonr, ttest_ind
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import brier_score_loss, roc_curve

# For survival
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

import warnings
warnings.filterwarnings('ignore')

def main():
    out_dir = 'outputs/comprehensive_analysis'
    os.makedirs(out_dir, exist_ok=True)
    
    # ---------------------------------------------------------
    # 0. Data Loading & Definition
    # ---------------------------------------------------------
    print("[1/5] Loading Data & Defining X, y...")
    pheno_df = pd.read_csv('BRCA/BRCA_phenotype.txt', sep='\t', index_col='idx')
    rna_df = pd.read_csv('BRCA/BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt', sep='\t', index_col=0).T
    surv_df = pd.read_csv('BRCA/BRCA_survival.txt', sep='\t', index_col=0)
    
    pheno_df.index = pheno_df.index.astype(str)
    rna_df.index = rna_df.index.astype(str)
    surv_df.index = surv_df.index.astype(str)
    
    median_immune = pheno_df['ESTIMATE_ImmuneScore'].median()
    pheno_df['Target_y'] = (pheno_df['ESTIMATE_ImmuneScore'] > median_immune).astype(int)
    
    df_merged = rna_df.join(pheno_df['Target_y'], how='inner').join(surv_df, how='inner').dropna(subset=['Target_y', 'OS_days', 'OS_event'])
    
    X_full = df_merged[rna_df.columns]
    y = df_merged['Target_y'].values
    
    # ---------------------------------------------------------
    # 1. Feature Selection (Correlation + LASSO)
    # ---------------------------------------------------------
    print("[2/5] Step 1: Feature Selection...")
    corrs = []
    for col in X_full.columns:
        if X_full[col].std() == 0: continue
        r, p = pearsonr(X_full[col], y)
        if p < 0.01:
            corrs.append((col, abs(r)))
            
    corrs.sort(key=lambda x: x[1], reverse=True)
    top_candidates = [x[0] for x in corrs[:300]]
    X_filtered = X_full[top_candidates]
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_filtered)
    
    lasso = LogisticRegressionCV(cv=5, penalty='l1', solver='liblinear', random_state=42)
    lasso.fit(X_scaled, y)
    
    coefs = lasso.coef_[0]
    selected_indices = np.where(coefs != 0)[0]
    
    if len(selected_indices) < 3 or len(selected_indices) > 20:
        top_indices = np.argsort(np.abs(coefs))[-10:]
        selected_genes = X_filtered.columns[top_indices].tolist()
        valid_coefs = coefs[top_indices]
    else:
        selected_genes = X_filtered.columns[selected_indices].tolist()
        valid_coefs = coefs[selected_indices]
        
    print(f"Selected {len(selected_genes)} Biomarkers: {selected_genes}")
    
    plt.figure(figsize=(10, 6))
    plt.barh(selected_genes, valid_coefs, color='skyblue')
    plt.title('LASSO Selected Biomarkers (Target: High Immune Subtype)')
    plt.xlabel('L1 Coefficient')
    plt.tight_layout()
    plt.savefig(f'{out_dir}/1_lasso_coefficients.png', dpi=300)
    plt.close()
    
    # Correlation Heatmap of selected features + Target
    X_sel = df_merged[selected_genes]
    corr_df = X_sel.copy()
    corr_df['Target'] = y
    corr_matrix = corr_df.corr(method='spearman')
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(corr_matrix, annot=True, fmt=".2f", cmap='coolwarm', center=0)
    plt.title("Correlation: Selected Biomarkers vs Target")
    plt.tight_layout()
    plt.savefig(f'{out_dir}/1_correlation_heatmap.png', dpi=300)
    plt.close()
    
    # ---------------------------------------------------------
    # 2. Diagnostic Model Construction
    # ---------------------------------------------------------
    print("[3/5] Step 2: Diagnostic Modeling...")
    lr_model = LogisticRegression(random_state=42)
    lr_model.fit(scaler.fit_transform(X_sel), y)
    model_scores = lr_model.predict_proba(scaler.transform(X_sel))[:, 1]
    df_merged['Model_Score'] = model_scores
    df_merged['Predicted_Group'] = (model_scores > 0.5).astype(int)
    
    dt_model = DecisionTreeClassifier(max_depth=3, random_state=42)
    dt_model.fit(X_sel, y)
    dot_data = export_graphviz(dt_model, out_file=None, feature_names=selected_genes, 
                               class_names=['Low Immune', 'High Immune'], filled=True, rounded=True)
    try:
        resp = requests.post("https://quickchart.io/graphviz", json={"graph": dot_data})
        if resp.status_code == 200:
            with open(f'{out_dir}/2_decision_tree.svg', 'wb') as f:
                f.write(resp.content)
    except:
        pass

    # ---------------------------------------------------------
    # 3. Validation of Biomarkers
    # ---------------------------------------------------------
    print("[4/5] Step 3: Biomarker Validation...")
    
    # a. Clustering
    sns.clustermap(X_sel.T, cmap='viridis', figsize=(10, 8), standard_scale=0)
    plt.title("Hierarchical Clustering of Selected Biomarkers")
    plt.savefig(f'{out_dir}/3a_clustering.png', dpi=300)
    plt.close()
    
    # b. t-test of Model Score by group
    plt.figure(figsize=(6, 4))
    sns.boxplot(x=y, y=model_scores, palette='Set2')
    t_stat, p_val = ttest_ind(model_scores[y==0], model_scores[y==1])
    plt.title(f"Model Score by Group (t-test p={p_val:.2e})")
    plt.xticks([0, 1], ['Low Immune', 'High Immune'])
    plt.ylabel('Predicted Probability')
    plt.tight_layout()
    plt.savefig(f'{out_dir}/3b_model_score_ttest.png', dpi=300)
    plt.close()
    
    # c. t-test of individual biomarkers
    fig, axes = plt.subplots(2, 5, figsize=(20, 8))
    axes = axes.flatten()
    for i, gene in enumerate(selected_genes):
        if i >= 10: break
        sns.boxplot(x=y, y=X_sel[gene], ax=axes[i], palette='Set2')
        _, p = ttest_ind(X_sel[gene][y==0], X_sel[gene][y==1])
        axes[i].set_title(f"{gene}\np={p:.2e}")
        axes[i].set_xticks([0, 1])
        axes[i].set_xticklabels(['Low', 'High'])
    plt.tight_layout()
    plt.savefig(f'{out_dir}/3c_biomarkers_ttest.png', dpi=300)
    plt.close()
    
    # d. Survival Analysis based on Predicted Group
    T = df_merged['OS_days']
    E = df_merged['OS_event']
    plt.figure(figsize=(8, 6))
    kmf = KaplanMeierFitter()
    if sum(df_merged['Predicted_Group'] == 1) > 0 and sum(df_merged['Predicted_Group'] == 0) > 0:
        kmf.fit(T[df_merged['Predicted_Group'] == 1], E[df_merged['Predicted_Group'] == 1], label='Predicted High Immune')
        kmf.plot(color='red')
        kmf.fit(T[df_merged['Predicted_Group'] == 0], E[df_merged['Predicted_Group'] == 0], label='Predicted Low Immune')
        kmf.plot(color='blue')
        res = logrank_test(T[df_merged['Predicted_Group'] == 1], T[df_merged['Predicted_Group'] == 0], 
                           event_observed_A=E[df_merged['Predicted_Group'] == 1], event_observed_B=E[df_merged['Predicted_Group'] == 0])
        plt.text(0.05, 0.1, f"Log-rank p={res.p_value:.4f}", transform=plt.gca().transAxes)
    plt.title("Survival Analysis by Predicted Biomarker Group")
    plt.savefig(f'{out_dir}/3d_survival_km.png', dpi=300)
    plt.close()
    
    # e. PCA & t-SNE of ONLY selected biomarkers
    pca = PCA(n_components=2).fit_transform(scaler.fit_transform(X_sel))
    tsne = TSNE(n_components=2, perplexity=10).fit_transform(scaler.fit_transform(X_sel))
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sns.scatterplot(x=pca[:,0], y=pca[:,1], hue=y, ax=axes[0], palette='Set1')
    axes[0].set_title('PCA of Selected Biomarkers')
    sns.scatterplot(x=tsne[:,0], y=tsne[:,1], hue=y, ax=axes[1], palette='Set1')
    axes[1].set_title('t-SNE of Selected Biomarkers')
    plt.tight_layout()
    plt.savefig(f'{out_dir}/3e_pca_tsne.png', dpi=300)
    plt.close()
    
    # f. DCA & Calibration
    plt.figure(figsize=(14, 6))
    plt.subplot(1, 2, 1)
    from sklearn.calibration import calibration_curve
    prob_true, prob_pred = calibration_curve(y, model_scores, n_bins=5)
    plt.plot(prob_pred, prob_true, marker='o')
    plt.plot([0, 1], [0, 1], linestyle='--')
    plt.title('Calibration Curve')
    
    plt.subplot(1, 2, 2)
    thresholds = np.linspace(0.01, 0.99, 100)
    net_benefits = []
    for thresh in thresholds:
        tp = np.sum((model_scores >= thresh) & (y == 1))
        fp = np.sum((model_scores >= thresh) & (y == 0))
        n = len(y)
        nb = (tp / n) - (fp / n) * (thresh / (1 - thresh))
        net_benefits.append(nb)
    plt.plot(thresholds, net_benefits, label='Model')
    plt.plot(thresholds, np.maximum(0, np.mean(y) - (1 - np.mean(y)) * (thresholds / (1 - thresholds))), label='Treat All')
    plt.plot(thresholds, np.zeros_like(thresholds), label='Treat None')
    plt.title('Decision Curve Analysis (DCA)')
    plt.ylim(-0.05, 0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{out_dir}/3f_dca_calibration.png', dpi=300)
    plt.close()
    
    # g. Volcano Plot
    print("[5/5] Generating Volcano Plot...")
    log_fc = []
    p_vals = []
    for col in X_full.columns:
        if X_full[col].std() == 0:
            log_fc.append(0)
            p_vals.append(1)
            continue
        fc = X_full[col][y==1].mean() - X_full[col][y==0].mean()
        _, p = ttest_ind(X_full[col][y==0], X_full[col][y==1])
        log_fc.append(fc)
        p_vals.append(p)
        
    p_vals = np.array(p_vals)
    p_vals[p_vals == 0] = 1e-300
    neg_log_p = -np.log10(p_vals)
    
    plt.figure(figsize=(10, 8))
    plt.scatter(log_fc, neg_log_p, color='gray', alpha=0.5)
    
    # Highlight selected genes
    sel_idx = [X_full.columns.get_loc(g) for g in selected_genes]
    plt.scatter(np.array(log_fc)[sel_idx], np.array(neg_log_p)[sel_idx], color='red', s=100, label='Selected Biomarkers')
    
    plt.title('Volcano Plot highlighting Selected Biomarkers')
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 P-value')
    plt.legend()
    plt.axhline(y=-np.log10(0.05), color='k', linestyle='--')
    plt.savefig(f'{out_dir}/3g_volcano_plot.png', dpi=300)
    plt.close()

    print(f"✅ Single Stream Pipeline Completed! Check '{out_dir}/' for all results.")

if __name__ == "__main__":
    main()
