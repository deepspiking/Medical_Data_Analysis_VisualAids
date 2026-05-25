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
    # 1. Feature Selection Phase 1 (Correlation + LASSO)
    # ---------------------------------------------------------
    print("[2/5] Step 1: Initial Feature Selection (LASSO)...")
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
        lasso_genes = X_filtered.columns[top_indices].tolist()
    else:
        lasso_genes = X_filtered.columns[selected_indices].tolist()
        
    # ---------------------------------------------------------
    # 2. Diagnostic Model Construction & Final Biomarker Extraction
    # ---------------------------------------------------------
    print("[3/5] Step 2: Diagnostic Modeling & Final Biomarker Extraction...")
    X_lasso = df_merged[lasso_genes]
    
    # Train Decision Tree on LASSO genes
    dt_model = DecisionTreeClassifier(max_depth=3, random_state=42)
    dt_model.fit(X_lasso, y)
    
    # Extract ONLY the genes actually used in the Decision Tree splits
    tree = dt_model.tree_
    used_indices = set(tree.feature[tree.feature != -2]) # -2 is _tree.TREE_UNDEFINED
    final_biomarkers = [lasso_genes[i] for i in used_indices]
    
    print(f"Final Biomarkers (Used in DT): {final_biomarkers}")
    X_final = df_merged[final_biomarkers]
    X_final_scaled = StandardScaler().fit_transform(X_final)
    
    # Train Logistic Regression ONLY on Final Biomarkers
    lr_model = LogisticRegression(random_state=42)
    lr_model.fit(X_final_scaled, y)
    
    # Save LR Coefficients Table
    coef_df = pd.DataFrame({
        'Biomarker': final_biomarkers,
        'LR_Coefficient': lr_model.coef_[0],
        'Odds_Ratio': np.exp(lr_model.coef_[0])
    }).sort_values(by='LR_Coefficient', ascending=False)
    coef_df.to_csv(f'{out_dir}/2_lr_coefficients.csv', index=False)
    
    # Create Table Image
    fig, ax = plt.subplots(figsize=(6, 2 + len(final_biomarkers)*0.5))
    ax.axis('tight')
    ax.axis('off')
    table_data = coef_df.round(4).values.tolist()
    table = ax.table(cellText=table_data, colLabels=coef_df.columns, cellLoc='center', loc='center')
    table.scale(1, 2)
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    plt.title('Logistic Regression Coefficients of Final Biomarkers', pad=20, fontsize=14, fontweight='bold')
    plt.savefig(f'{out_dir}/2_lr_coefficients_table.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Generate Model Scores
    lr_scores = lr_model.predict_proba(X_final_scaled)[:, 1]
    dt_scores = dt_model.predict_proba(X_lasso)[:, 1]  # DT trained on X_lasso, but only uses final_biomarkers internally
    
    df_merged['LR_Score'] = lr_scores
    df_merged['DT_Score'] = dt_scores
    df_merged['LR_Predicted_Group'] = (lr_scores > 0.5).astype(int)
    df_merged['DT_Predicted_Group'] = dt_model.predict(X_lasso)
    
    # Save DT Image
    dot_data = export_graphviz(dt_model, out_file=None, feature_names=lasso_genes, 
                               class_names=['Low Immune', 'High Immune'], filled=True, rounded=True)
    try:
        resp = requests.post("https://quickchart.io/graphviz", json={"graph": dot_data})
        if resp.status_code == 200:
            with open(f'{out_dir}/2_decision_tree.svg', 'wb') as f:
                f.write(resp.content)
    except:
        pass

    # ---------------------------------------------------------
    # 3. Validation of FINAL Biomarkers
    # ---------------------------------------------------------
    print("[4/5] Step 3: Validation of FINAL Biomarkers...")
    
    # a. Correlation Heatmap of Final Biomarkers
    corr_df = X_final.copy()
    corr_df['Target'] = y
    corr_matrix = corr_df.corr(method='spearman')
    plt.figure(figsize=(8, 6))
    sns.heatmap(corr_matrix, annot=True, fmt=".2f", cmap='coolwarm', center=0)
    plt.title("Correlation: Final Biomarkers vs Target")
    plt.tight_layout()
    plt.savefig(f'{out_dir}/3a_correlation_heatmap.png', dpi=300)
    plt.close()

    # b. Clustering of Final Biomarkers
    sns.clustermap(X_final.T, cmap='viridis', figsize=(8, 6), standard_scale=0)
    plt.title("Hierarchical Clustering of Final Biomarkers", pad=20)
    plt.savefig(f'{out_dir}/3b_clustering.png', dpi=300)
    plt.close()
    
    # c. t-test of Model Score by group (BOTH LR and DT)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    sns.boxplot(x=y, y=lr_scores, palette='Set2', ax=axes[0])
    t_stat_lr, p_val_lr = ttest_ind(lr_scores[y==0], lr_scores[y==1])
    axes[0].set_title(f"Logistic Regression Score (p={p_val_lr:.2e})")
    axes[0].set_xticks([0, 1])
    axes[0].set_xticklabels(['Low Immune', 'High Immune'])
    axes[0].set_ylabel('Predicted Probability')
    
    sns.boxplot(x=y, y=dt_scores, palette='Set2', ax=axes[1])
    t_stat_dt, p_val_dt = ttest_ind(dt_scores[y==0], dt_scores[y==1])
    axes[1].set_title(f"Decision Tree Leaf Score (p={p_val_dt:.2e})")
    axes[1].set_xticks([0, 1])
    axes[1].set_xticklabels(['Low Immune', 'High Immune'])
    axes[1].set_ylabel('Leaf Probability')
    
    plt.tight_layout()
    plt.savefig(f'{out_dir}/3c_model_scores_ttest.png', dpi=300)
    plt.close()
    
    # d. t-test of individual final biomarkers
    n_cols = len(final_biomarkers)
    fig, axes = plt.subplots(1, n_cols, figsize=(n_cols*4, 5))
    if n_cols == 1: axes = [axes]
    for i, gene in enumerate(final_biomarkers):
        sns.boxplot(x=y, y=X_final[gene], ax=axes[i], palette='Set2')
        _, p = ttest_ind(X_final[gene][y==0], X_final[gene][y==1])
        axes[i].set_title(f"{gene}\np={p:.2e}")
        axes[i].set_xticks([0, 1])
        axes[i].set_xticklabels(['Low', 'High'])
    plt.tight_layout()
    plt.savefig(f'{out_dir}/3d_biomarkers_ttest.png', dpi=300)
    plt.close()
    
    # e. Survival Analysis based on Predicted Group (Using LR Predictions)
    T = df_merged['OS_days']
    E = df_merged['OS_event']
    plt.figure(figsize=(8, 6))
    kmf = KaplanMeierFitter()
    if sum(df_merged['LR_Predicted_Group'] == 1) > 0 and sum(df_merged['LR_Predicted_Group'] == 0) > 0:
        kmf.fit(T[df_merged['LR_Predicted_Group'] == 1], E[df_merged['LR_Predicted_Group'] == 1], label='LR High Immune')
        kmf.plot(color='red')
        kmf.fit(T[df_merged['LR_Predicted_Group'] == 0], E[df_merged['LR_Predicted_Group'] == 0], label='LR Low Immune')
        kmf.plot(color='blue')
        res = logrank_test(T[df_merged['LR_Predicted_Group'] == 1], T[df_merged['LR_Predicted_Group'] == 0], 
                           event_observed_A=E[df_merged['LR_Predicted_Group'] == 1], event_observed_B=E[df_merged['LR_Predicted_Group'] == 0])
        plt.text(0.05, 0.1, f"Log-rank p={res.p_value:.4f}", transform=plt.gca().transAxes)
    plt.title("Survival Analysis by LR Predicted Group")
    plt.savefig(f'{out_dir}/3e_survival_km.png', dpi=300)
    plt.close()
    
    # f. PCA & t-SNE of ONLY FINAL biomarkers
    if len(final_biomarkers) > 1:
        pca = PCA(n_components=2).fit_transform(X_final_scaled)
        tsne = TSNE(n_components=2, perplexity=10, random_state=42).fit_transform(X_final_scaled)
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        sns.scatterplot(x=pca[:,0], y=pca[:,1], hue=y, ax=axes[0], palette='Set1', s=80)
        axes[0].set_title('PCA of Final Biomarkers')
        sns.scatterplot(x=tsne[:,0], y=tsne[:,1], hue=y, ax=axes[1], palette='Set1', s=80)
        axes[1].set_title('t-SNE of Final Biomarkers')
        plt.tight_layout()
        plt.savefig(f'{out_dir}/3f_pca_tsne.png', dpi=300)
        plt.close()
    
    # g. DCA & Calibration for BOTH models
    plt.figure(figsize=(14, 6))
    
    # Calibration Curve
    plt.subplot(1, 2, 1)
    from sklearn.calibration import calibration_curve
    prob_true_lr, prob_pred_lr = calibration_curve(y, lr_scores, n_bins=5)
    prob_true_dt, prob_pred_dt = calibration_curve(y, dt_scores, n_bins=5)
    plt.plot(prob_pred_lr, prob_true_lr, marker='o', label='Logistic Regression')
    plt.plot(prob_pred_dt, prob_true_dt, marker='s', label='Decision Tree')
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray')
    plt.title('Calibration Curve')
    plt.legend()
    
    # DCA Plot
    plt.subplot(1, 2, 2)
    thresholds = np.linspace(0.01, 0.99, 100)
    def calc_nb(scores, y_true, thresh):
        tp = np.sum((scores >= thresh) & (y_true == 1))
        fp = np.sum((scores >= thresh) & (y_true == 0))
        n = len(y_true)
        return (tp / n) - (fp / n) * (thresh / (1 - thresh))
        
    nb_lr = [calc_nb(lr_scores, y, t) for t in thresholds]
    nb_dt = [calc_nb(dt_scores, y, t) for t in thresholds]
    nb_all = np.maximum(0, np.mean(y) - (1 - np.mean(y)) * (thresholds / (1 - thresholds)))
    
    plt.plot(thresholds, nb_lr, label='Logistic Regression')
    plt.plot(thresholds, nb_dt, label='Decision Tree')
    plt.plot(thresholds, nb_all, label='Treat All', linestyle=':')
    plt.plot(thresholds, np.zeros_like(thresholds), label='Treat None', color='k')
    plt.title('Decision Curve Analysis (DCA)')
    plt.ylim(-0.05, 0.6)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(f'{out_dir}/3g_dca_calibration.png', dpi=300)
    plt.close()
    
    # h. Volcano Plot highlighting FINAL biomarkers
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
    plt.scatter(log_fc, neg_log_p, color='lightgray', alpha=0.5)
    
    sel_idx = [X_full.columns.get_loc(g) for g in final_biomarkers]
    plt.scatter(np.array(log_fc)[sel_idx], np.array(neg_log_p)[sel_idx], color='red', s=150, edgecolor='black', label='Final Biomarkers (DT+LR)')
    
    # Add text labels to highlighted points
    for i, g in zip(sel_idx, final_biomarkers):
        plt.text(log_fc[i]+0.1, neg_log_p[i]+0.1, g, fontsize=10, weight='bold')
    
    plt.title('Volcano Plot highlighting Final Biomarkers')
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 P-value')
    plt.legend()
    plt.axhline(y=-np.log10(0.05), color='k', linestyle='--')
    plt.savefig(f'{out_dir}/3h_volcano_plot.png', dpi=300)
    plt.close()

    print(f"✅ Single Stream Pipeline Completed! Check '{out_dir}/' for all results.")

if __name__ == "__main__":
    main()
