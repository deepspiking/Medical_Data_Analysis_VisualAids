import os
import subprocess
import time

def print_header(title):
    print("\n" + "="*60)
    print(f"🚀 {title}")
    print("="*60)

def run_script(script_name, step_num, total_steps):
    print(f"\n[Step {step_num}/{total_steps}] Executing {script_name}...")
    start_time = time.time()
    
    script_path = os.path.join("src", script_name)
    if not os.path.exists(script_path):
        print(f"❌ Error: {script_path} not found!")
        return False
        
    try:
        result = subprocess.run(["python", script_path], capture_output=True, text=True)
        if result.returncode == 0:
            elapsed = time.time() - start_time
            print(f"✅ Success! (Took {elapsed:.2f}s)")
            return True
        else:
            print(f"❌ Failed with error code {result.returncode}")
            print(result.stderr)
            return False
    except Exception as e:
        print(f"❌ Execution failed: {e}")
        return False

def main():
    print_header("Medical Data Analysis & Visualization Pipeline")
    print("Starting the end-to-end orchestration of all 11 analytical modules.")
    
    scripts = [
        "01_comparative_analysis.py",
        "02_categorical_correlation.py",
        "03_survival_analysis.py",
        "04_predictive_modeling.py",
        "05_concordance_analysis.py",
        "06_causal_inference.py",
        "07_advanced_visualization.py",
        "08_hierarchical_clustering.py",
        "09_dimensionality_reduction.py",
        "10_pathway_enrichment.py",
        "11_feature_selection_lasso.py"
    ]
    
    total = len(scripts)
    success_count = 0
    
    for i, script in enumerate(scripts, 1):
        if run_script(script, i, total):
            success_count += 1
            
    print_header("Pipeline Execution Summary")
    print(f"Total Steps: {total}")
    print(f"Successful:  {success_count}")
    print(f"Failed:      {total - success_count}")
    
    if success_count == total:
        print("\n🎉 All analyses completed successfully! Please check the 'outputs/' directory for all plots.")
    else:
        print("\n⚠️ Pipeline finished with some errors. Please check the logs above.")

if __name__ == "__main__":
    main()
