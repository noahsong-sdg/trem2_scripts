#!/usr/bin/env python3

import os
import subprocess
import sys
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# --- Configuration ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RECEPTOR_FILE = os.path.join(SCRIPT_DIR, "../data/receptor/cluster1_fixed.pdb")
LIGAND_DIR = os.path.join(SCRIPT_DIR, "../data/column_one/ligands_sdf_split/")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "../results/c1_outputs/")

# Box parameters (same as original script)
CENTER_X, CENTER_Y, CENTER_Z = 42.328, 28.604, 21.648
SIZE_X, SIZE_Y, SIZE_Z = 22.5, 22.5, 22.5

# Chunk processing parameters (for cluster time limit resilience)
LIGANDS_PER_CHUNK = 1000  # Number of ligands per mcdock execution

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
    "--savedir": os.path.join(OUTPUT_DIR, "multiconfdockresult"),
    "--batch_size": "500", # 1200 caused broken pipe and ran out of memory
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


def validate_sdf_file(sdf_path):
    """
    Validate an SDF file using Open Babel to check if it's readable.
    
    Args:
        sdf_path (str): Path to the SDF file to validate
        
    Returns:
        bool: True if file is valid, False otherwise
    """
    try:
        # Use Open Babel to validate the SDF file
        result = subprocess.run(
            ['obabel', sdf_path, '-osdf', '-O', '/dev/null'],
            capture_output=True,
            text=True,
            timeout=30
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
        return False

def filter_valid_ligands(ligand_files):
    """
    Filter out invalid/corrupted SDF files before processing.
    
    Args:
        ligand_files (list): List of ligand file paths
        
    Returns:
        tuple: (valid_files, invalid_files)
    """
    valid_files = []
    invalid_files = []
    
    logging.info(f"Validating {len(ligand_files)} ligand files...")
    
    for i, ligand_file in enumerate(ligand_files, 1):
        if i % 100 == 0:
            logging.info(f"Validated {i}/{len(ligand_files)} files...")
        
        if validate_sdf_file(ligand_file):
            valid_files.append(ligand_file)
        else:
            invalid_files.append(ligand_file)
            logging.warning(f"Invalid SDF file detected: {os.path.basename(ligand_file)}")
    
    logging.info(f"Validation complete: {len(valid_files)} valid, {len(invalid_files)} invalid")
    return valid_files, invalid_files

def get_completed_ligands(output_dir):
    """
    Check which ligands already have output files in multiconfdockresult directory.
    
    Args:
        output_dir (str): Base output directory containing multiconfdockresult subdirectory
        
    Returns:
        set: Set of ligand names that already have SDF output files
    """
    result_dir = os.path.join(output_dir, "multiconfdockresult")
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


def create_chunks(ligand_files, chunk_size):
    """
    Split ligand files into chunks of specified size.
    
    Args:
        ligand_files (list): List of ligand file paths
        chunk_size (int): Number of ligands per chunk
        
    Returns:
        list: List of chunks, where each chunk is a list of ligand file paths
    """
    chunks = []
    for i in range(0, len(ligand_files), chunk_size):
        chunk = ligand_files[i:i + chunk_size]
        chunks.append(chunk)
    return chunks


def run_mcdock_chunk(chunk_ligands, chunk_num, total_chunks):
    """
    Run UniDock mcdock on a single chunk of ligands.
    
    Args:
        chunk_ligands (list): List of ligand file paths for this chunk
        chunk_num (int): Current chunk number (1-indexed)
        total_chunks (int): Total number of chunks
        
    Returns:
        bool: True if chunk completed successfully, False otherwise
    """
    logging.info(f"\n=== Processing Chunk {chunk_num}/{total_chunks} ({len(chunk_ligands)} ligands) ===")
    logging.info(f"UniDock will process these in batches of 800 ligands internally")
    
    # Create chunk-specific ligand list file
    chunk_ligand_list = os.path.join(OUTPUT_DIR, f"chunk_{chunk_num}_ligand_list.txt")
    with open(chunk_ligand_list, "w") as f:
        for ligand in chunk_ligands:
            f.write(f"{os.path.abspath(ligand)}\n")
    
    # Build command for this chunk
    cmd = ["unidocktools", "mcdock"]
    for flag, value in MCDOCK_FLAGS.items():
        if flag == "--gen_conf":
            if value is not None:
                cmd.append(flag)
        elif value is not None:
            cmd.extend([flag, value])
    cmd.extend(["--ligand_index", chunk_ligand_list])
    
    logging.info(f"Running chunk {chunk_num} command:")
    
    logging.info(" ".join(cmd))
    
    try:
        result = subprocess.run(cmd, check=True, text=True, capture_output=True)
        logging.info(f"Chunk {chunk_num} completed successfully")
        logging.info(result.stdout)
        
        # Clean up chunk ligand list file
        os.remove(chunk_ligand_list)
        return True
        
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in chunk {chunk_num}: {e}")
        if e.stdout:
            logging.error(f"[stdout] {e.stdout}")
        if e.stderr:
            logging.error(f"[stderr] {e.stderr}")
        return False
    except KeyboardInterrupt:
        logging.error(f"\nChunk {chunk_num} interrupted by user!")
        raise


def reset_progress():
    """Reset progress by removing the multiconfdockresult directory."""
    result_dir = os.path.join(OUTPUT_DIR, "multiconfdockresult")
    if os.path.exists(result_dir):
        import shutil
        shutil.rmtree(result_dir)
        logging.info(f"Reset complete: Removed {result_dir}")
        return True
    else:
        logging.info(f"No existing results found at {result_dir}")
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
        logging.error(f"Error: Receptor file not found at {RECEPTOR_FILE}")
        exit(1)
    if not os.path.exists(LIGAND_DIR):
        logging.error(f"Error: Ligand directory not found at {LIGAND_DIR}")
        exit(1)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Collect all ligand files (SDF or PDBQT) recursively
    all_ligand_files = sorted(
        [str(p) for p in Path(LIGAND_DIR).rglob("*.sdf")] +
        [str(p) for p in Path(LIGAND_DIR).rglob("*.pdbqt")]
    )
    if not all_ligand_files:
        logging.error(f"No ligand files (.sdf or .pdbqt) found in {LIGAND_DIR} or its subdirectories")
        exit(1)

    # Check for completed ligands and filter
    logging.info(f"Found {len(all_ligand_files)} total ligand files")
    completed_ligands = get_completed_ligands(OUTPUT_DIR)
    
    if completed_ligands:
        logging.info(f"Resume detected: {len(completed_ligands)} ligands already completed")
        ligand_files = filter_completed_ligands(all_ligand_files, completed_ligands)
        logging.info(f"Remaining ligands to process: {len(ligand_files)}")
        
        if not ligand_files:
            logging.info("All ligands have already been processed!")
            logging.info("Use --reset flag to start over from scratch")
            exit(0)
    else:
        logging.info("Starting fresh docking run")
        ligand_files = all_ligand_files

    # Validate ligand files to filter out corrupted ones
    valid_ligands, invalid_ligands = filter_valid_ligands(ligand_files)
    
    if invalid_ligands:
        logging.warning(f"Found {len(invalid_ligands)} invalid SDF files that will be skipped:")
        for invalid_file in invalid_ligands[:10]:  # Show first 10
            logging.warning(f"  - {os.path.basename(invalid_file)}")
        if len(invalid_ligands) > 10:
            logging.warning(f"  ... and {len(invalid_ligands) - 10} more")
    
    if not valid_ligands:
        logging.error("No valid ligand files found after validation!")
        exit(1)
    
    ligand_files = valid_ligands  # Use only valid files
    logging.info(f"Proceeding with {len(ligand_files)} valid ligand files")

    # Create chunks for processing
    chunks = create_chunks(ligand_files, LIGANDS_PER_CHUNK)
    total_chunks = len(chunks)
    
    logging.info(f"\nChunk Processing Plan:")
    logging.info(f"  Total ligands to process: {len(ligand_files)}")
    logging.info(f"  Ligands per chunk: {LIGANDS_PER_CHUNK}")
    logging.info(f"  Total chunks: {total_chunks}")
    logging.info(f"  UniDock batch size: 800 ligands (internal)")
    
    # Process chunks sequentially
    successful_chunks = 0
    failed_chunks = 0
    
    try:
        for chunk_num, chunk_ligands in enumerate(chunks, 1):
            # Check if any ligands in this chunk are already completed
            chunk_completed = get_completed_ligands(OUTPUT_DIR)
            remaining_in_chunk = filter_completed_ligands(chunk_ligands, chunk_completed)
            
            if not remaining_in_chunk:
                logging.info(f"\nSkipping Chunk {chunk_num}/{total_chunks} - all ligands already completed")
                successful_chunks += 1
                continue
            
            if len(remaining_in_chunk) < len(chunk_ligands):
                logging.info(f"Chunk {chunk_num}: {len(remaining_in_chunk)}/{len(chunk_ligands)} ligands remaining")
            
            # Run the chunk with remaining ligands
            chunk_success = run_mcdock_chunk(remaining_in_chunk, chunk_num, total_chunks)
            
            if chunk_success:
                successful_chunks += 1
                
                # Show progress after each chunk
                current_completed = get_completed_ligands(OUTPUT_DIR)
                total_completed = len(current_completed)
                logging.info(f"Progress: {total_completed}/{len(all_ligand_files)} ligands completed")
            else:
                failed_chunks += 1
                logging.error(f"Chunk {chunk_num} failed - stopping execution")
                break
                
    except KeyboardInterrupt:
        logging.error(f"\n\nProcessing interrupted by user!")
        current_completed = get_completed_ligands(OUTPUT_DIR)
        if current_completed:
            logging.info(f"Progress saved: {len(current_completed)} ligands completed")
            logging.info("Run the script again to resume from where it left off")
        exit(1)
    
    # Final summary
    final_completed = get_completed_ligands(OUTPUT_DIR)
    total_completed = len(final_completed)
    total_ligands = len(all_ligand_files)
    
    logging.info(f"\n=== FINAL SUMMARY ===")
    logging.info(f"Chunks processed: {successful_chunks}/{total_chunks}")
    logging.info(f"Total ligands completed: {total_completed}/{total_ligands}")
    logging.info(f"Output directory: {os.path.join(OUTPUT_DIR, 'multiconfdockresult')}")
    
    if failed_chunks > 0:
        logging.error(f"Failed chunks: {failed_chunks}")
        logging.info("Run the script again to retry failed chunks")
        exit(1)
    elif total_completed == total_ligands:
        logging.info("All ligands completed successfully!")
    else:
        remaining = total_ligands - total_completed
        logging.info(f"Remaining ligands: {remaining}")
        logging.info("Run the script again to resume from where it left off")

if __name__ == "__main__":
    main() 
