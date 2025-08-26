import os
import pandas as pd
import glob
import re

# --- Configuration ---
DOCKING_OUTPUT_DIR = "../results/c1_outputs/mcresult/"
ANALYSIS_RESULTS_FILE = "../results/mcdock_c1p1_summary.csv"

def parse_unidock_sdf(filepath):
    """
    Parses a UniDock output SDF file to extract docking scores.
    
    Args:
        filepath (str): Path to the Uni-Dock output SDF file.
    
    Returns:
        list: List of dictionaries with conformation info and scores.
    """
    results = []
    try:
        with open(filepath, 'r') as f:
            content = f.read()
            
        # Split by $$$$ to get individual conformations
        conformations = content.split('$$$$')
        
        for i, conf in enumerate(conformations):
            if not conf.strip():
                continue
                
            # Extract ligand name from first line
            lines = conf.strip().split('\n')
            if not lines:
                continue
                
            ligand_name = lines[0].strip()
            
            # Find Uni-Dock RESULT section
            result_match = re.search(r'ENERGY=\s*([-\d.]+)', conf)
            if result_match:
                energy = float(result_match.group(1))
                
                # Extract bounds if available
                lower_bound_match = re.search(r'LOWER_BOUND=\s*([-\d.]+)', conf)
                upper_bound_match = re.search(r'UPPER_BOUND=\s*([-\d.]+)', conf)
                
                lower_bound = float(lower_bound_match.group(1)) if lower_bound_match else None
                upper_bound = float(upper_bound_match.group(1)) if upper_bound_match else None
                
                results.append({
                    'ligand_name': ligand_name,
                    'conformation': i + 1,
                    'energy': energy,
                    'lower_bound': lower_bound,
                    'upper_bound': upper_bound,
                    'file_path': filepath
                })
                
    except Exception as e:
        print(f"Error parsing file {filepath}: {e}")
        return []
    
    return results

def main():
    print("--- Analyzing UniDock Results ---")
    
    all_results = []
    
    # Look for SDF files in the output directory
    search_pattern = os.path.join(DOCKING_OUTPUT_DIR, "*.sdf")
    sdf_files = glob.glob(search_pattern, recursive=True)
    
    if not sdf_files:
        print(f"No SDF output files found in {DOCKING_OUTPUT_DIR}")
        print("Please check the DOCKING_OUTPUT_DIR path.")
        exit()
    
    print(f"Found {len(sdf_files)} SDF files to analyze.")
    
    for sdf_file in sdf_files:
        print(f"Processing: {sdf_file}")
        file_results = parse_unidock_sdf(sdf_file)
        all_results.extend(file_results)
    
    if not all_results:
        print("No results were processed. Exiting.")
        exit()
    
    # Create DataFrame
    results_df = pd.DataFrame(all_results)
    
    # Sort by energy (most negative is best)
    results_df = results_df.sort_values(by='energy', ascending=True)
    
    # Save results
    os.makedirs(os.path.dirname(ANALYSIS_RESULTS_FILE), exist_ok=True)
    results_df.to_csv(ANALYSIS_RESULTS_FILE, index=False)
    
    print(f"\n--- Analysis Complete ---")
    print(f"Summary saved to: {ANALYSIS_RESULTS_FILE}")
    print(f"Total conformations analyzed: {len(results_df)}")
    
    # Show top results
    print("\nTop 1000 Docking Results:")
    print(results_df.head(1000)[['ligand_name', 'conformation', 'energy']].to_string())
    
    # Summary statistics
    print(f"\nSummary Statistics:")
    print(f"Best energy: {results_df['energy'].min():.3f}")
    print(f"Average energy: {results_df['energy'].mean():.3f}")
    print(f"Number of unique ligands: {results_df['ligand_name'].nunique()}")

if __name__ == "__main__":
    main()
