import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import os

def main():
    os.makedirs('outputs/survival_plots', exist_ok=True)
    
    survival_file = "BRCA/BRCA_survival.txt"
    risk_file = "outputs/pipeline_state/patient_risk.csv"
    
    if not os.path.exists(survival_file) or not os.path.exists(risk_file):
        print("Required files not found. Ensure pipeline script 04 ran first.")
        return

    print("Loading survival and ML risk data...")
    df_surv = pd.read_csv(survival_file, sep='\t', index_col=0)
    df_risk = pd.read_csv(risk_file, index_col=0)
    
    df = df_surv.join(df_risk, how='inner').dropna(subset=['OS_days', 'OS_event'])
    
    print("1. Kaplan-Meier Survival Analysis (High vs Low Risk from ML Model)")
    T = df['OS_days']
    E = df['OS_event']
    
    high_mask = df['Risk_Group'] == 'High Risk'
    low_mask = df['Risk_Group'] == 'Low Risk'
    
    plt.figure(figsize=(8, 6))
    kmf_high = KaplanMeierFitter()
    if high_mask.sum() > 0:
        kmf_high.fit(T[high_mask], event_observed=E[high_mask], label='High Risk')
        kmf_high.plot_survival_function(color='red')
        
    kmf_low = KaplanMeierFitter()
    if low_mask.sum() > 0:
        kmf_low.fit(T[low_mask], event_observed=E[low_mask], label='Low Risk')
        kmf_low.plot_survival_function(color='blue')
        
    plt.title('Kaplan-Meier Survival Curve by ML-Predicted Risk')
    plt.xlabel('Days')
    plt.ylabel('Survival Probability')
    
    # Log-rank test
    if high_mask.sum() > 0 and low_mask.sum() > 0:
        results = logrank_test(T[high_mask], T[low_mask], event_observed_A=E[high_mask], event_observed_B=E[low_mask])
        p_val = results.p_value
        plt.text(0.05, 0.2, f"Log-rank p={p_val:.4f}", transform=plt.gca().transAxes, fontsize=12)
        print(f"Log-rank p-value: {p_val:.4f}")
        
    plt.savefig('outputs/survival_plots/km_survival_all.png', dpi=300)
    plt.close()
    print("Saved 'outputs/survival_plots/km_survival_all.png'.")

if __name__ == "__main__":
    main()
