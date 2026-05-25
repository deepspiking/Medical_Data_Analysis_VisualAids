import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import numpy as np

# Patch for older sklearn versions
np.float = float
np.int = int
np.bool = bool
np.object = object

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier, plot_tree
from sklearn.metrics import accuracy_score, classification_report
from sklearn.preprocessing import StandardScaler

def load_data(pheno_path, surv_path):
    # Load phenotype data (features)
    pheno_df = pd.read_csv(pheno_path, sep='\t', index_col='idx')
    
    # Load survival data (outcome)
    surv_df = pd.read_csv(surv_path, sep='\t', index_col='case_id')
    
    # Merge on patient ID
    df = pheno_df.join(surv_df, how='inner')
    
    # Use CIBERSORT and PROGENy features as a proxy for RNAseq/proteomics derived numerical features
    feature_cols = [col for col in df.columns if col.startswith('CIBERSORT_') or col.startswith('PROGENy_')]
    
    # Target: Overall Survival event
    # Drop rows with missing target
    df = df.dropna(subset=['OS_event'])
    
    X = df[feature_cols].copy()
    y = df['OS_event'].astype(int)
    
    # Handle any missing values in features by filling with median
    X = X.fillna(X.median())
    
    return X, y

def multimodal_fusion_template(X_rna, X_prot):
    """
    A dummy mock template for multimodal fusion using deep learning.
    This illustrates how one might combine RNAseq and Proteomics features
    using a neural network in a PyTorch-like pseudo-code structure.
    """
    print("\n--- Multimodal Fusion Deep Learning Template ---")
    template_code = '''
import torch
import torch.nn as nn

class MultimodalFusionNet(nn.Module):
    def __init__(self, rna_input_dim, prot_input_dim, hidden_dim=64):
        super(MultimodalFusionNet, self).__init__()
        
        # RNA branch
        self.rna_net = nn.Sequential(
            nn.Linear(rna_input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.3)
        )
        
        # Proteomics branch
        self.prot_net = nn.Sequential(
            nn.Linear(prot_input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.3)
        )
        
        # Fusion layer
        self.fusion_layer = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1),
            nn.Sigmoid()
        )
        
    def forward(self, rna_x, prot_x):
        rna_feat = self.rna_net(rna_x)
        prot_feat = self.prot_net(prot_x)
        
        # Concatenate features
        fused = torch.cat((rna_feat, prot_feat), dim=1)
        
        out = self.fusion_layer(fused)
        return out
        
# Example usage:
# model = MultimodalFusionNet(rna_input_dim=100, prot_input_dim=100)
# predictions = model(rna_tensor, prot_tensor)
'''
    print(template_code)

def main():
    pheno_path = 'BRCA/BRCA_phenotype.txt'
    surv_path = 'BRCA/BRCA_survival.txt'
    
    if not os.path.exists(pheno_path) or not os.path.exists(surv_path):
        print(f"Data files not found. Ensure {pheno_path} and {surv_path} exist.")
        return
        
    print("Loading data...")
    X, y = load_data(pheno_path, surv_path)
    
    print(f"Dataset shape: {X.shape}, Target shape: {y.shape}")
    print(f"Class distribution:\n{y.value_counts()}")
    
    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Standardize features for Logistic Regression
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # 1. Logistic Regression
    print("\n--- Training Logistic Regression ---")
    lr_model = LogisticRegression(max_iter=1000, random_state=42)
    lr_model.fit(X_train_scaled, y_train)
    
    lr_preds = lr_model.predict(X_test_scaled)
    print("Logistic Regression Accuracy:", accuracy_score(y_test, lr_preds))
    print(classification_report(y_test, lr_preds, zero_division=0))
    
    # 2. Decision Tree
    print("\n--- Training Decision Tree ---")
    dt_model = DecisionTreeClassifier(max_depth=4, random_state=42, class_weight='balanced')
    dt_model.fit(X_train, y_train)
    
    dt_preds = dt_model.predict(X_test)
    print("Decision Tree Accuracy:", accuracy_score(y_test, dt_preds))
    print(classification_report(y_test, dt_preds, zero_division=0))
    
    # Save Decision Tree Visualization using Graphviz API for a beautiful layout
    from sklearn.tree import export_graphviz
    import requests
    
    output_dir = 'outputs/predictive_modeling'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    dt_viz_path = os.path.join(output_dir, 'decision_tree_viz.svg')
    
    short_feature_names = [col.replace('CIBERSORT_', '').replace('PROGENy_', '') for col in X.columns]
    
    dot_data = export_graphviz(dt_model, out_file=None, 
                               feature_names=short_feature_names, 
                               class_names=['Alive', 'Deceased'], 
                               filled=True, rounded=True, special_characters=True)
                               
    try:
        # Use an external API to render the graphviz dot string into a beautiful PNG image
        # This completely avoids the overlapping issues of matplotlib's plot_tree
        resp = requests.post("https://quickchart.io/graphviz", json={"graph": dot_data}, timeout=10)
        if resp.status_code == 200:
            with open(dt_viz_path, 'wb') as f:
                f.write(resp.content)
            print(f"\nDecision tree visualization (Beautiful Graphviz) saved to {dt_viz_path}")
        else:
            raise Exception(f"API Error {resp.status_code}")
    except Exception as e:
        print(f"Graphviz API rendering failed ({e}), falling back to matplotlib...")
        fig, ax = plt.subplots(figsize=(12, 8))
        plot_tree(dt_model, feature_names=short_feature_names, class_names=['Alive', 'Deceased'], 
                  filled=True, rounded=True, ax=ax)
        plt.title("Decision Tree for BRCA Survival Prediction")
        plt.savefig(dt_viz_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"\nDecision tree visualization (Fallback) saved to {dt_viz_path}")
    
    # 3. Multimodal Fusion Template
    multimodal_fusion_template(None, None)

if __name__ == "__main__":
    main()
