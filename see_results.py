import os
import pandas as pd
import glob

# --- Configuration ---
DOCKING_OUTPUT_DIR = "../results/outputs/"
ANALYSIS_RESULTS_FILE = "../results/summary.csv"

# --- Helper Function to Parse Uni-Dock Output Files ---
def parse_unidock_pdbqt_vina(filepath):
    """
    Parses a Uni-Dock output PDBQT file (Vina-like format) to extract docking scores.
    Assumes the scores are in lines starting with "REMARK VINA RESULT:".

    Args:
        filepath (str): Path to the Uni-Dock output PDBQT file for a single ligand.

    Returns:
        list: A list of dictionaries, where each dictionary contains 
              {'mode': mode_number, 'affinity_kcal_mol': affinity} for each binding mode.
              Returns an empty list if no scores are found or the file cannot be parsed.
    """
    scores = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.split()
                    # Example line: REMARK VINA RESULT:      -7.5      0.000      0.000
                    # We need the affinity (2nd field after REMARK VINA RESULT:)
                    # and the mode number (implicitly the order they appear)
                    if len(parts) >= 4:
                        try:
                            affinity = float(parts[3])
                            # Mode number is 1-based from the order of appearance
                            scores.append({"affinity_kcal_mol": affinity})
                        except ValueError:
                            print(f"Warning: Could not parse score from line in {filepath}: {line.strip()}")
                elif line.startswith("MODEL") and scores: # Vina PDBQT separates modes by MODEL records
                    # If we hit a new MODEL and have scores, assign mode numbers to previous set
                    for i, score_entry in enumerate(scores):
                        if "mode" not in score_entry: # Assign mode number if not already set
                             score_entry["mode"] = i + 1 
                    
    except FileNotFoundError:
        print(f"Error: File not found {filepath}")
        return []
    except Exception as e:
        print(f"Error parsing file {filepath}: {e}")
        return []
    
    # Final assignment of mode numbers if not done by MODEL keyword (e.g. if only one mode)
    for i, score_entry in enumerate(scores):
        if "mode" not in score_entry:
            score_entry["mode"] = i + 1

    return scores


if __name__ == "__main__":
    print("--- Analyzing Uni-Dock Results ---")

    all_results = []

    search_pattern = os.path.join(DOCKING_OUTPUT_DIR, "*_docked.pdbqt")

    docked_files = glob.glob(search_pattern, recursive=True)

    if not docked_files:
        print(f"No PDBQT output files found in {DOCKING_OUTPUT_DIR} matching pattern {search_pattern}.")
        print("Please check the DOCKING_OUTPUT_DIR and the Uni-Dock execution.")
        exit()

    print(f"Found {len(docked_files)} PDBQT files to analyze.")

    for pdbqt_file in docked_files:
        ligand_name_with_ext = os.path.basename(pdbqt_file)
        ligand_name = os.path.splitext(ligand_name_with_ext)[0].replace("_docked", "") # Basic name extraction
        
        # If files are in subdirectories named after the ligand:
        # ligand_name = os.path.basename(os.path.dirname(pdbqt_file))

        print(f"Processing: {pdbqt_file} (Ligand: {ligand_name})")
        
        # Parse the PDBQT file for scores
        # This function needs to be adapted based on the exact format of Uni-Dock's output PDBQT
        # For Vina-like output, scores are in REMARK VINA RESULT lines.
        ligand_scores = parse_unidock_pdbqt_vina(pdbqt_file)

        if not ligand_scores:
            print(f"Warning: No scores found or error parsing for {ligand_name} in {pdbqt_file}")
            all_results.append({
                "ligand_name": ligand_name,
                "file_path": pdbqt_file,
                "mode": None,
                "affinity_kcal_mol": None,
                "error": "Parsing failed or no scores"
            })
            continue

        for score_info in ligand_scores:
            all_results.append({
                "ligand_name": ligand_name,
                "file_path": pdbqt_file,
                "mode": score_info.get("mode"),
                "affinity_kcal_mol": score_info.get("affinity_kcal_mol")
            })

    if not all_results:
        print("No results were processed. Exiting.")
        exit()

    # Create a Pandas DataFrame for easy analysis and saving
    results_df = pd.DataFrame(all_results)

    # Sort by affinity (most negative is best)
    results_df = results_df.sort_values(by="affinity_kcal_mol", ascending=True)

    # Save the summary to a CSV file
    os.makedirs(os.path.dirname(ANALYSIS_RESULTS_FILE), exist_ok=True)
    results_df.to_csv(ANALYSIS_RESULTS_FILE, index=False)

    print(f"\n--- Analysis Complete ---")
    print(f"Summary of results saved to: {ANALYSIS_RESULTS_FILE}")

    # Display top N results (e.g., top 10)
    print("\nTop 10 Docking Results (Best Affinity per Mode):")
    # If multiple modes per ligand, you might want to see the best mode per ligand first.
    best_mode_per_ligand = results_df.loc[results_df.groupby("ligand_name")["affinity_kcal_mol"].idxmin()]
    print(best_mode_per_ligand.head(10).to_string())

    print("\nFull sorted results are in the CSV file.")
