import os
import pandas as pd
import glob
import re

# --- Configuration ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DOCKING_OUTPUT_DIR = os.path.join(SCRIPT_DIR, "../results/c1_outputs/mcresult")
ANALYSIS_RESULTS_FILE = os.path.join(SCRIPT_DIR, "../results/mcdock_c1p1_summary.csv")

def parse_docking_sdf(filepath):
    """
    Parses a docking output SDF file to extract docking scores.
    
    Args:
        filepath (str): Path to the docking output SDF file.
    
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
            
            # Look for docking_score field
            score_match = re.search(r'>\s*<docking_score>\s*\((\d+)\)\s*\n([-\d.]+)', conf)
            if score_match:
                conformation_num = int(score_match.group(1))
                score = float(score_match.group(2))
                
                results.append({
                    'ligand_name': ligand_name,
                    'conformation': conformation_num,
                    'docking_score': score,
                    'file_path': filepath
                })
                
    except Exception as e:
        print(f"Error parsing file {filepath}: {e}")
        return []
    
    return results

def main():
    print("--- Analyzing Docking Results ---")
    print(f"Searching in: {DOCKING_OUTPUT_DIR}")
    
    all_results = []
    
    # Look for SDF files in the output directory
    search_pattern = os.path.join(DOCKING_OUTPUT_DIR, "*.sdf")
    sdf_files = glob.glob(search_pattern, recursive=True)
    
    if not sdf_files:
        print(f"No SDF output files found in {DOCKING_OUTPUT_DIR}")
        print("Please check the DOCKING_OUTPUT_DIR path.")
        exit()
    
    print(f"Found {len(sdf_files)} SDF files to analyze.")
    
    # Process first few files to test
    test_files = sdf_files[:5]
    print(f"Testing with first {len(test_files)} files...")
    
    for sdf_file in test_files:
        print(f"Processing: {os.path.basename(sdf_file)}")
        file_results = parse_docking_sdf(sdf_file)
        if file_results:
            print(f"  Found {len(file_results)} conformations")
            all_results.extend(file_results)
        else:
            print(f"  No results found")
    
    if not all_results:
        print("No results were processed from test files. Exiting.")
        print("Let's examine the first file structure:")
        
        if sdf_files:
            first_file = sdf_files[0]
            print(f"\nExamining: {first_file}")
            try:
                with open(first_file, 'r') as f:
                    content = f.read()
                    # Show the first 500 characters
                    print("First 500 characters:")
                    print(content[:500])
                    print("\n" + "="*50)
                    # Look for docking_score patterns
                    score_matches = re.findall(r'>\s*<docking_score>.*?\n([-\d.]+)', content, re.DOTALL)
                    print(f"Found {len(score_matches)} docking_score entries")
                    if score_matches:
                        print("Sample scores:", score_matches[:5])
            except Exception as e:
                print(f"Error reading file: {e}")
        exit()
    
    # If test was successful, process all files
    print(f"\nTest successful! Processing all {len(sdf_files)} files...")
    
    for i, sdf_file in enumerate(sdf_files):
        if i % 1000 == 0:
            print(f"Progress: {i}/{len(sdf_files)} files processed")
        file_results = parse_docking_sdf(sdf_file)
        all_results.extend(file_results)
    
    if not all_results:
        print("No results were processed. Exiting.")
        exit()
    
    # Create DataFrame
    results_df = pd.DataFrame(all_results)
    
    # Sort by docking_score (most negative is best)
    results_df = results_df.sort_values(by='docking_score', ascending=True)
    
    # Save results
    os.makedirs(os.path.dirname(ANALYSIS_RESULTS_FILE), exist_ok=True)
    results_df.to_csv(ANALYSIS_RESULTS_FILE, index=False)
    
    print(f"\n--- Analysis Complete ---")
    print(f"Summary saved to: {ANALYSIS_RESULTS_FILE}")
    print(f"Total conformations analyzed: {len(results_df)}")
    
    # Show top results
    print("\nTop 20 Docking Results:")
    print(results_df.head(20)[['ligand_name', 'conformation', 'docking_score']].to_string())
    
    # Summary statistics
    print(f"\nSummary Statistics:")
    print(f"Best docking score: {results_df['docking_score'].min():.3f}")
    print(f"Average docking score: {results_df['docking_score'].mean():.3f}")
    print(f"Number of unique ligands: {results_df['ligand_name'].nunique()}")

if __name__ == "__main__":
    main()
