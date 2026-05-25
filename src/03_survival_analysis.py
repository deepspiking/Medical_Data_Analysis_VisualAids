import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import os

def main():
    # File paths
    survival_file = "BRCA/BRCA_survival.txt"
    phenotype_file = "BRCA/BRCA_phenotype.txt"

    # Ensure files exist
    if not os.path.exists(survival_file) or not os.path.exists(phenotype_file):
        print("Data files not found. Ensure you are running this from the correct directory.")
        return

    print("Loading data...")
    # Load data
    df_surv = pd.read_csv(survival_file, sep='\t')
    df_pheno = pd.read_csv(phenotype_file, sep='\t')
    
    # Rename phenotype idx to case_id for merging
    df_pheno = df_pheno.rename(columns={'idx': 'case_id'})
    
    print("Merging data...")
    # Merge on case_id
    df = pd.merge(df_surv, df_pheno, on='case_id', how='inner')
    
    # We will use Overall Survival (OS_days, OS_event)
    df = df.dropna(subset=['OS_days', 'OS_event'])
    T = df['OS_days']
    E = df['OS_event']
    
    print("1. Kaplan-Meier Survival Analysis")
    kmf = KaplanMeierFitter()
    kmf.fit(T, event_observed=E, label='All Patients')
    
    plt.figure(figsize=(8, 6))
    kmf.plot_survival_function()
    plt.title('Kaplan-Meier Survival Curve (Overall Survival)')
    plt.xlabel('Days')
    plt.ylabel('Survival Probability')
    os.makedirs('outputs/survival_plots', exist_ok=True)
    plt.savefig('outputs/survival_plots/km_survival_all.png')
    plt.close()
    print("Saved 'outputs/survival_plots/km_survival_all.png'.")
    
    print("\n2. Log-Rank Test")
    # Compare groups based on the median of ESTIMATE_ImmuneScore
    score_col = 'ESTIMATE_ImmuneScore'
    if score_col in df.columns:
        median_score = df[score_col].median()
        
        # Split into high and low groups
        high_mask = df[score_col] > median_score
        low_mask = ~high_mask
        
        T_high = df[high_mask]['OS_days']
        E_high = df[high_mask]['OS_event']
        
        T_low = df[low_mask]['OS_days']
        E_low = df[low_mask]['OS_event']
        
        # Plot comparative KM curves
        plt.figure(figsize=(8, 6))
        kmf_high = KaplanMeierFitter()
        kmf_high.fit(T_high, event_observed=E_high, label=f'High {score_col}')
        kmf_high.plot_survival_function()
        
        kmf_low = KaplanMeierFitter()
        kmf_low.fit(T_low, event_observed=E_low, label=f'Low {score_col}')
        kmf_low.plot_survival_function()
        
        plt.title(f'Survival by {score_col}')
        plt.xlabel('Days')
        plt.ylabel('Survival Probability')
        plt.savefig(f'outputs/survival_plots/km_survival_{score_col.lower()}.png')
        plt.close()
        print(f"Saved 'outputs/survival_plots/km_survival_{score_col.lower()}.png'.")
        
        # Log-rank test
        results = logrank_test(T_high, T_low, event_observed_A=E_high, event_observed_B=E_low)
        print(f"Log-rank test for {score_col} (High vs Low):")
        results.print_summary()
    else:
        print(f"Column '{score_col}' not found for Log-rank test.")
        
    print("\n3. Cox Proportional Hazards Model")
    # Select a few numeric covariates for the Cox model
    covariates = ['OS_days', 'OS_event', 'ESTIMATE_ImmuneScore', 'ESTIMATE_StromalScore', 'TMB']
    
    # Filter available covariates and drop NaNs
    available_covariates = [c for c in covariates if c in df.columns]
    cox_df = df[available_covariates].dropna()
    
    # Ensure variables are numeric
    for col in available_covariates:
        cox_df[col] = pd.to_numeric(cox_df[col], errors='coerce')
    cox_df = cox_df.dropna()
    
    if len(cox_df) > 10:  # arbitrary small threshold to ensure enough data
        cph = CoxPHFitter()
        # Fit the model
        cph.fit(cox_df, duration_col='OS_days', event_col='OS_event')
        print("Cox Proportional Hazards Model Summary:")
        cph.print_summary()
    else:
        print("Not enough data to fit the Cox PH model.")
        
if __name__ == "__main__":
    main()
