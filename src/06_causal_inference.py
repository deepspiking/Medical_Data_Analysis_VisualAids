import numpy as np

# Monkey-patch deprecated np.float to fix scikit-learn version mismatch
np.float = float
np.int = int
np.object = object
np.bool = bool

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors
import os

# Create results dir if needed
os.makedirs("results", exist_ok=True)

def generate_mock_clinical_data(n=1000):
    """Generates mock clinical data with confounding variables."""
    np.random.seed(42)
    
    # Covariates
    age = np.random.normal(55, 12, n)
    bmi = np.random.normal(28, 6, n)
    blood_pressure = np.random.normal(130, 20, n)
    
    # Treatment assignment (confounded by covariates)
    # Higher age and bmi increase likelihood of treatment
    logit_p = -3 + 0.04 * age + 0.08 * bmi - 0.01 * blood_pressure
    p_treatment = 1 / (1 + np.exp(-logit_p))
    treatment = np.random.binomial(1, p_treatment)
    
    # Outcome (affected by treatment and covariates)
    outcome = 20 + 5 * treatment + 0.3 * age + 0.4 * bmi + np.random.normal(0, 3, n)
    
    return pd.DataFrame({
        'Age': age,
        'BMI': bmi,
        'Blood_Pressure': blood_pressure,
        'Treatment': treatment,
        'Outcome': outcome
    })

def calc_smd(df, covariates, treatment_col, weights=None):
    """Calculates Absolute Standardized Mean Differences (ASMD)."""
    smds = []
    treated = df[df[treatment_col] == 1]
    control = df[df[treatment_col] == 0]
    
    for cov in covariates:
        if weights is None:
            mean_t = treated[cov].mean()
            mean_c = control[cov].mean()
            var_t = treated[cov].var()
            var_c = control[cov].var()
        else:
            w_t = df.loc[df[treatment_col] == 1, weights]
            w_c = df.loc[df[treatment_col] == 0, weights]
            
            mean_t = np.average(treated[cov], weights=w_t)
            mean_c = np.average(control[cov], weights=w_c)
            
            var_t = np.average((treated[cov] - mean_t)**2, weights=w_t)
            var_c = np.average((control[cov] - mean_c)**2, weights=w_c)
            
        pooled_std = np.sqrt((var_t + var_c) / 2)
        smd = np.abs(mean_t - mean_c) / pooled_std
        smds.append(smd)
        
    return smds

def main():
    print("Generating mock clinical data...")
    data = generate_mock_clinical_data(2500)
    covariates = ['Age', 'BMI', 'Blood_Pressure']
    
    print(f"Dataset generated: {len(data)} patients ({data['Treatment'].sum()} treated, {len(data)-data['Treatment'].sum()} control)")
    
    # 1. Calculate Propensity Scores
    print("Calculating Propensity Scores...")
    X = data[covariates]
    y = data['Treatment']
    
    model = LogisticRegression(random_state=42)
    model.fit(X, y)
    data['Propensity_Score'] = model.predict_proba(X)[:, 1]
    
    # 2. Propensity Score Matching (Nearest Neighbor, 1:1)
    print("Performing Propensity Score Matching...")
    treated = data[data['Treatment'] == 1].copy()
    control = data[data['Treatment'] == 0].copy()
    
    nn = NearestNeighbors(n_neighbors=1)
    nn.fit(control[['Propensity_Score']])
    distances, indices = nn.kneighbors(treated[['Propensity_Score']])
    
    matched_control = control.iloc[indices.flatten()]
    matched_data = pd.concat([treated, matched_control])
    
    # 3. Inverse Probability of Treatment Weighting (IPTW)
    print("Calculating IPTW weights...")
    data['IPTW'] = np.where(
        data['Treatment'] == 1,
        1 / data['Propensity_Score'],
        1 / (1 - data['Propensity_Score'])
    )
    
    # Trim extreme weights (1st and 99th percentiles) to avoid instability
    lower_q, upper_q = data['IPTW'].quantile([0.01, 0.99])
    data['IPTW_Trimmed'] = np.clip(data['IPTW'], lower_q, upper_q)
    
    # 4. Calculate Absolute Standardized Mean Differences (ASMD)
    print("Calculating Standardized Mean Differences for Love Plot...")
    
    unadjusted_smd = calc_smd(data, covariates, 'Treatment')
    matched_smd = calc_smd(matched_data, covariates, 'Treatment')
    iptw_smd = calc_smd(data, covariates, 'Treatment', weights='IPTW_Trimmed')
    
    # 5. Plot Love Plot
    print("Generating Love Plot...")
    smd_df = pd.DataFrame({
        'Covariate': covariates,
        'Unadjusted': unadjusted_smd,
        'PS Matched (1:1)': matched_smd,
        'IPTW (Trimmed)': iptw_smd
    })
    
    smd_df_melted = smd_df.melt(id_vars='Covariate', var_name='Method', value_name='ASMD')
    
    plt.figure(figsize=(10, 6))
    sns.set_theme(style="whitegrid")
    
    # Using specific markers for better readability
    markers = {"Unadjusted": "X", "PS Matched (1:1)": "o", "IPTW (Trimmed)": "s"}
    
    ax = sns.scatterplot(
        data=smd_df_melted, 
        y='Covariate', 
        x='ASMD', 
        hue='Method', 
        style='Method',
        markers=markers,
        s=120,
        alpha=0.8
    )
    
    plt.axvline(x=0.1, color='red', linestyle='--', linewidth=1.5, label='Threshold (0.1)')
    plt.title('Love Plot: Covariate Balance Before and After Adjustment', fontsize=14)
    plt.xlabel('Absolute Standardized Mean Difference (ASMD)', fontsize=12)
    plt.ylabel('')
    
    # Get handles and labels for legend, add threshold
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles=handles, labels=labels, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.grid(axis='x', linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    os.makedirs('outputs/causal_inference', exist_ok=True)
    output_path = 'outputs/causal_inference/love_plot.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Success: Love plot saved to {output_path}")

if __name__ == "__main__":
    main()
