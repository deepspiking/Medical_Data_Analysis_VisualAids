import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.metrics import cohen_kappa_score
import statsmodels.api as sm

def calculate_cohens_kappa(rater1, rater2):
    """
    Calculate Cohen's Kappa for categorical agreement between two raters.
    """
    kappa = cohen_kappa_score(rater1, rater2)
    print(f"Cohen's Kappa Score: {kappa:.4f}")
    return kappa

def plot_bland_altman(measure1, measure2, save_path=None):
    """
    Generate a Bland-Altman plot using statsmodels to assess agreement between two continuous measures.
    """
    f, ax = plt.subplots(1, figsize=(8, 5))
    sm.graphics.mean_diff_plot(measure1, measure2, ax=ax)
    plt.title("Bland-Altman Plot")
    
    if save_path:
        plt.savefig(save_path)
        print(f"Bland-Altman plot saved to {save_path}")
    else:
        print("Bland-Altman plot generated.")
        
    plt.close()

def main():
    print("--- Concordance Analysis ---")
    
    # 1. Cohen's Kappa Example (Categorical Data)
    print("\n1. Categorical Agreement (Cohen's Kappa)")
    # Mock data: 2 raters, classifying 50 samples into 3 categories (0, 1, 2)
    np.random.seed(42)
    rater_A = np.random.choice([0, 1, 2], size=50, p=[0.5, 0.3, 0.2])
    # Make rater B somewhat agree with rater A
    rater_B = rater_A.copy()
    noise_idx = np.random.choice(50, size=15, replace=False)
    rater_B[noise_idx] = np.random.choice([0, 1, 2], size=15)
    
    calculate_cohens_kappa(rater_A, rater_B)
    
    # 2. Bland-Altman Example (Continuous Data)
    print("\n2. Continuous Agreement (Bland-Altman Plot)")
    # Mock data: 2 devices measuring a continuous variable (e.g., tumor size)
    measure_device_1 = np.random.normal(loc=50, scale=10, size=100)
    # Device 2 has a slight bias and some noise compared to Device 1
    measure_device_2 = measure_device_1 + np.random.normal(loc=2, scale=3, size=100)
    
    import os
    os.makedirs('outputs/concordance_plots', exist_ok=True)
    plot_path = "outputs/concordance_plots/bland_altman_plot.png"
    plot_bland_altman(measure_device_1, measure_device_2, save_path=plot_path)

if __name__ == "__main__":
    main()
