#!/usr/bin/env python3
"""
Minimal UniDock docking script for all ligands in a directory.
- Centralizes all CLI flags for easy modification.
- Supports both SDF and PDBQT ligands (no SDF-only enforcement).
- No batching, no SLURM, no tranche logic.
"""
import os
import subprocess
from pathlib import Path

# --- Configuration ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RECEPTOR_FILE = os.path.join(SCRIPT_DIR, "../data/receptor/cluster1_fixed.pdb")
LIGAND_DIR = os.path.join(SCRIPT_DIR, "../data/column_one/ligands_sdf_split/")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "../results/mcdock_outputs_column_one/")

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


def main():
    # Validate input files/directories
    if not os.path.exists(RECEPTOR_FILE):
        print(f"Error: Receptor file not found at {RECEPTOR_FILE}")
        exit(1)
    if not os.path.exists(LIGAND_DIR):
        print(f"Error: Ligand directory not found at {LIGAND_DIR}")
        exit(1)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Collect all ligand files (SDF or PDBQT) recursively
    ligand_files = sorted(
        [str(p) for p in Path(LIGAND_DIR).rglob("*.sdf")] +
        [str(p) for p in Path(LIGAND_DIR).rglob("*.pdbqt")]
    )
    if not ligand_files:
        print(f"No ligand files (.sdf or .pdbqt) found in {LIGAND_DIR} or its subdirectories")
        exit(1)

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
    except subprocess.CalledProcessError as e:
        print(f"Error running UniDock mcdock: {e}")
        if e.stdout:
            print("\n[stdout]\n", e.stdout)
        if e.stderr:
            print("\n[stderr]\n", e.stderr)
        exit(1)

if __name__ == "__main__":
    main() 
