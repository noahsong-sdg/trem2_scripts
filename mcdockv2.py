#!/usr/bin/env python3

import os
import subprocess
import sys
from pathlib import Path

# --- Configuration ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RECEPTOR_FILE = os.path.join(SCRIPT_DIR, "../data/receptor/cluster1_fixed.pdb")
LIGAND_DIR = os.path.join(SCRIPT_DIR, "../data/column_one/ligands_sdf_split/")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "../results/c1_outputs/")

# Box parameters (same as original script)
CENTER_X, CENTER_Y, CENTER_Z = 42.328, 28.604, 21.648
SIZE_X, SIZE_Y, SIZE_Z = 22.5, 22.5, 22.5

# Centralized UniDock mcdock flags (add/modify here)
MCDOCK_FLAGS = {
    "--receptor": RECEPTOR_FILE,
    # Ligands will be set dynamically
    "--center_x": str(CENTER_X),
    "--center_y": str(CENTER_Y),
    "--center_z": str(CENTER_Z),
    "--size_x": str(SIZE_X),
    "--size_y": str(SIZE_Y),
    "--size_z": str(SIZE_Z),
    # "--workdir": os.path.join(OUTPUT_DIR, "MultiConfDock"),
    # "--savedir": os.path.join(OUTPUT_DIR, "MultiConfDock-Result"),
    "--batch_size": "800", # 1200 caused broken pipe and ran out of memory
    "--scoring_function_rigid_docking": "vina",
    "--exhaustiveness_rigid_docking": "32",
    "--num_modes_rigid_docking": "3",
    "--topn_rigid_docking": "20",
    "--scoring_function_local_refine": "vina",
    "--exhaustiveness_local_refine": "64",
    "--num_modes_local_refine": "1",
    "--topn_local_refine": "1",
    "--min_rmsd": "0.3",
    "--max_num_confs_per_ligand": "50",
    "--gen_conf": "1",  # Flag only, no value
}


def get_completed_ligands(output_dir):
    """
    Check which ligands already have output files in MultiConfDock-Result directory.
    
    Args:
        output_dir (str): Base output directory containing MultiConfDock-Result subdirectory
        
    Returns:
        set: Set of ligand names that already have SDF output files
    """
    result_dir = os.path.join(output_dir, "MultiConfDock-Result")
    if not os.path.exists(result_dir):
        return set()
    
    completed = set()
    for sdf_file in os.listdir(result_dir):
        if sdf_file.endswith('.sdf') and os.path.getsize(os.path.join(result_dir, sdf_file)) > 0:
            # Extract ligand name (remove .sdf extension)
            ligand_name = sdf_file[:-4]
            completed.add(ligand_name)
    return completed


def filter_completed_ligands(ligand_files, completed_ligands):
    """
    Remove ligands that are already completed from the input list.
    
    Args:
        ligand_files (list): List of ligand file paths
        completed_ligands (set): Set of completed ligand names
        
    Returns:
        list: Filtered list of ligand files that still need processing
    """
    remaining = []
    for ligand_file in ligand_files:
        ligand_name = os.path.splitext(os.path.basename(ligand_file))[0]
        if ligand_name not in completed_ligands:
            remaining.append(ligand_file)
    return remaining


def reset_progress():
    """Reset progress by removing the MultiConfDock-Result directory."""
    result_dir = os.path.join(OUTPUT_DIR, "MultiConfDock-Result")
    if os.path.exists(result_dir):
        import shutil
        shutil.rmtree(result_dir)
        print(f"Reset complete: Removed {result_dir}")
        return True
    else:
        print(f"No existing results found at {result_dir}")
        return False


def main():
    # Check for command line arguments
    if len(sys.argv) > 1:
        if sys.argv[1] == "--reset":
            reset_progress()
            exit(0)
        elif sys.argv[1] == "--help" or sys.argv[1] == "-h":
            print("UniDock mcdock with resume capability")
            print("Usage:")
            print("  python mcdockv2.py          # Run/resume docking")
            print("  python mcdockv2.py --reset  # Reset progress and start over")
            exit(0)

    # Validate input files/directories
    if not os.path.exists(RECEPTOR_FILE):
        print(f"Error: Receptor file not found at {RECEPTOR_FILE}")
        exit(1)
    if not os.path.exists(LIGAND_DIR):
        print(f"Error: Ligand directory not found at {LIGAND_DIR}")
        exit(1)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Collect all ligand files (SDF or PDBQT) recursively
    all_ligand_files = sorted(
        [str(p) for p in Path(LIGAND_DIR).rglob("*.sdf")] +
        [str(p) for p in Path(LIGAND_DIR).rglob("*.pdbqt")]
    )
    if not all_ligand_files:
        print(f"No ligand files (.sdf or .pdbqt) found in {LIGAND_DIR} or its subdirectories")
        exit(1)

    # Check for completed ligands and filter
    print(f"Found {len(all_ligand_files)} total ligand files")
    completed_ligands = get_completed_ligands(OUTPUT_DIR)
    
    if completed_ligands:
        print(f"Resume detected: {len(completed_ligands)} ligands already completed")
        ligand_files = filter_completed_ligands(all_ligand_files, completed_ligands)
        print(f"Remaining ligands to process: {len(ligand_files)}")
        
        if not ligand_files:
            print("All ligands have already been processed!")
            print("Use --reset flag to start over from scratch")
            exit(0)
    else:
        print("Starting fresh docking run")
        ligand_files = all_ligand_files

    # Write ligand list to file (avoids long command lines)
    ligand_list_file = os.path.join(OUTPUT_DIR, "ligand_list.txt")
    with open(ligand_list_file, "w") as f:
        for ligand in ligand_files:
            f.write(f"{os.path.abspath(ligand)}\n")

    # Build command
    cmd = ["unidocktools", "mcdock"]
    for flag, value in MCDOCK_FLAGS.items():
        if flag == "--gen_conf":
            if value is not None:
                cmd.append(flag)
        elif value is not None:
            cmd.extend([flag, value])
    cmd.extend(["--ligand_index", ligand_list_file])

    print("Running UniDock mcdock with command:")
    print(" ", " ".join(cmd))

    try:
        result = subprocess.run(cmd, check=True, text=True, capture_output=True)
        print("\n=== UniDock mcdock output ===")
        print(result.stdout)
        print("\n=== UniDock mcdock completed successfully ===")
        
        # Final summary
        final_completed = get_completed_ligands(OUTPUT_DIR)
        total_completed = len(final_completed)
        total_ligands = len(all_ligand_files)
        
        print(f"\n=== FINAL SUMMARY ===")
        print(f"Total ligands processed: {total_completed}/{total_ligands}")
        print(f"Output directory: {os.path.join(OUTPUT_DIR, 'MultiConfDock-Result')}")
        
        if total_completed == total_ligands:
            print("All ligands completed successfully!")
        else:
            remaining = total_ligands - total_completed
            print(f"Remaining ligands: {remaining}")
            print("Run the script again to resume from where it left off")
            
    except subprocess.CalledProcessError as e:
        print(f"Error running UniDock mcdock: {e}")
        if e.stdout:
            print("\n[stdout]\n", e.stdout)
        if e.stderr:
            print("\n[stderr]\n", e.stderr)
        
        # Show resume information even on failure
        final_completed = get_completed_ligands(OUTPUT_DIR)
        if final_completed:
            print(f"\nProgress saved: {len(final_completed)} ligands completed")
            print("Run the script again to resume from where it left off")
        exit(1)
    except KeyboardInterrupt:
        print(f"\n\nInterrupted by user!")
        final_completed = get_completed_ligands(OUTPUT_DIR)
        if final_completed:
            print(f"Progress saved: {len(final_completed)} ligands completed")
            print("Run the script again to resume from where it left off")
        exit(1)

if __name__ == "__main__":
    main() 
