#!/usr/bin/env python3
"""
Performs a single UniDock docking calculation for one ligand.
Receives all paths and parameters as command-line arguments.
This script is designed to be called by a Slurm job array.
"""
import os
import subprocess
import argparse

def run_single_unidock(args):
    """Executes a single UniDock docking run."""

    # Ensure the output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Define the final output path for the docked ligand
    ligand_basename = os.path.splitext(os.path.basename(args.ligand_path))[0]
    output_path = os.path.join(args.output_dir, f"{ligand_basename}_out.pdbqt")

    # --- This is our new, simpler 'resume' functionality ---
    # If the output file already exists and is not empty, we assume it was
    # completed successfully in a previous run.
    if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
        print(f"INFO: Output {output_path} already exists. Skipping.")
        return 0 # Success code for skipping

    command = [
        "unidock",  # Assumes 'unidock' is in the system PATH
        "--receptor", args.receptor_path,
        "--ligand", args.ligand_path,
        "--dir", args.output_dir,
        "--center_x", str(args.center_x),
        "--center_y", str(args.center_y),
        "--center_z", str(args.center_z),
        "--size_x", str(args.size_x),
        "--size_y", str(args.size_y),
        "--size_z", str(args.size_z),
        "--scoring", "vinardo",
        "--num_modes", "5"
    ]

    print(f"Executing command: {' '.join(command)}")
    try:
        # We use subprocess.run to execute the command.
        # We capture output to prevent it from cluttering the main Slurm log unless there's an error.
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        print(f"Successfully docked {args.ligand_path}")
        return 0 # Success
    except subprocess.CalledProcessError as e:
        print(f"ERROR: UniDock failed for ligand: {args.ligand_path}")
        print(f"Return code: {e.returncode}")
        print(f"Stderr: {e.stderr}")
        return 1 # Failure

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run a single UniDock job.")
    # Input files
    parser.add_argument("--receptor_path", type=str, required=True, help="Path to the receptor PDBQT file.")
    parser.add_argument("--ligand_path", type=str, required=True, help="Path to the single ligand PDBQT file.")
    # Output directory
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save docking results.")
    # Binding site definition
    parser.add_argument("--center_x", type=float, required=True)
    parser.add_argument("--center_y", type=float, required=True)
    parser.add_argument("--center_z", type=float, required=True)
    parser.add_argument("--size_x", type=float, required=True)
    parser.add_argument("--size_y", type=float, required=True)
    parser.add_argument("--size_z", type=float, required=True)

    args = parser.parse_args()
    
    exit_code = run_single_unidock(args)
    exit(exit_code)
