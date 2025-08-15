#!/usr/bin/env python3

import os
import subprocess
import sys
import signal
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

# Global shutdown flag
shutdown_requested = False

def signal_handler(signum, frame):
    global shutdown_requested
    logging.warning(f"Signal {signum} received. Shutting down gracefully...")
    shutdown_requested = True

# Configuration
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RECEPTOR_FILE = os.path.join(SCRIPT_DIR, "../data/receptor/cluster1_fixed.pdb")
LIGAND_DIR = os.path.join(SCRIPT_DIR, "../data/column_two/ligands_sdf_split/")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "../results/c2_outputs/")

# Box parameters
CENTER_X, CENTER_Y, CENTER_Z = 42.328, 28.604, 21.648
SIZE_X, SIZE_Y, SIZE_Z = 22.5, 22.5, 22.5

# Processing parameters
LIGANDS_PER_CHUNK = 1000

# UniDock flags
MCDOCK_FLAGS = {
    "--receptor": RECEPTOR_FILE,
    "--center_x": str(CENTER_X), "--center_y": str(CENTER_Y), "--center_z": str(CENTER_Z),
    "--size_x": str(SIZE_X), "--size_y": str(SIZE_Y), "--size_z": str(SIZE_Z),
    "--savedir": os.path.join(OUTPUT_DIR, "mcresult"),
    "--batch_size": "250",
    "--scoring_function_rigid_docking": "vina",
    "--exhaustiveness_rigid_docking": "32",
    "--num_modes_rigid_docking": "3",
    "--topn_rigid_docking": "20",
    "--scoring_function_local_refine": "vina",
    "--exhaustiveness_local_refine": "32",
    "--num_modes_local_refine": "1",
    "--topn_local_refine": "1",
    "--min_rmsd": "0.3",
    "--max_num_confs_per_ligand": "20",
    "--gen_conf": "1",
}

def get_completed_ligands():
    """Get set of completed ligand names."""
    result_dir = os.path.join(OUTPUT_DIR, "mcresult")
    if not os.path.exists(result_dir):
        return set()
    
    completed = set()
    for sdf_file in os.listdir(result_dir):
        if sdf_file.endswith('.sdf') and os.path.getsize(os.path.join(result_dir, sdf_file)) > 0:
            completed.add(sdf_file[:-4])
    return completed

def filter_completed_ligands(ligand_files, completed_ligands):
    """Remove completed ligands from list."""
    return [f for f in ligand_files if os.path.splitext(os.path.basename(f))[0] not in completed_ligands]

def create_chunks(ligand_files, chunk_size):
    """Split files into chunks."""
    return [ligand_files[i:i + chunk_size] for i in range(0, len(ligand_files), chunk_size)]

def quick_file_test(ligand_files):
    """Quick test if files are accessible."""
    accessible = []
    for ligand_file in ligand_files:
        try:
            if os.path.exists(ligand_file) and os.path.getsize(ligand_file) > 0:
                with open(ligand_file, 'r', encoding='utf-8', errors='ignore') as f:
                    if f.readline().strip():
                        accessible.append(ligand_file)
        except:
            pass
    return accessible

def run_mcdock_chunk(chunk_ligands, chunk_num, total_chunks):
    """Run UniDock on a chunk of ligands."""
    logging.info(f"Chunk {chunk_num}/{total_chunks}: {len(chunk_ligands)} ligands")
    
    # Quick file test
    accessible_files = quick_file_test(chunk_ligands)
    if len(accessible_files) != len(chunk_ligands):
        logging.error(f"Chunk {chunk_num}: {len(chunk_ligands) - len(accessible_files)} inaccessible files")
        return False
    
    # Create ligand list file
    chunk_list = os.path.join(OUTPUT_DIR, f"chunk_{chunk_num}_ligand_list.txt")
    with open(chunk_list, "w") as f:
        for ligand in chunk_ligands:
            f.write(f"{os.path.abspath(ligand)}\n")
    
    # Build command
    cmd = ["unidocktools", "mcdock"]
    for flag, value in MCDOCK_FLAGS.items():
        if flag == "--gen_conf":
            cmd.append(flag)
        elif value is not None:
            cmd.extend([flag, value])
    cmd.extend(["--ligand_index", chunk_list])
    
    try:
        result = subprocess.run(cmd, check=True, text=True, capture_output=True, timeout=36000,
                               stdin=subprocess.DEVNULL)
        logging.info(f"Chunk {chunk_num} completed")
        os.remove(chunk_list)
        return True
        
    except subprocess.CalledProcessError as e:
        logging.error(f"Chunk {chunk_num} failed: {e}")
        if e.stderr:
            logging.error(f"STDERR: {e.stderr}")
        if e.stdout:
            logging.error(f"STDOUT: {e.stdout}")
        
        if "Bad input file" in str(e.stderr):
            logging.error("Bad input file error detected")
            # Try to extract the problematic filename
            error_lines = str(e.stderr).split('\n')
            for line in error_lines:
                if "Bad input file" in line and ".sdf" in line:
                    logging.error(f"Problematic file: {line}")
                    # Extract just the filename
                    if "obabel_" in line:
                        parts = line.split('/')
                        for part in parts:
                            if part.endswith('.sdf'):
                                logging.error(f"Failing SDF file: {part}")
                                break
        return False
    except subprocess.TimeoutExpired:
        logging.error(f"Chunk {chunk_num} timed out")
        return False

def reset_progress():
    """Reset progress by removing mcresult directory."""
    result_dir = os.path.join(OUTPUT_DIR, "mcresult")
    if os.path.exists(result_dir):
        import shutil
        shutil.rmtree(result_dir)
        logging.info("Progress reset")
        return True
    return False

def run_diagnostic_tests():
    """Quick diagnostic tests."""
    logging.info("=== DIAGNOSTIC TESTS ===")
    
    # Test 1: Check directories
    if not os.path.exists(RECEPTOR_FILE):
        logging.error(f"Receptor not found: {RECEPTOR_FILE}")
        return False
    if not os.path.exists(LIGAND_DIR):
        logging.error(f"Ligand dir not found: {LIGAND_DIR}")
        return False
    
    # Test 2: Find files
    all_files = sorted([str(p) for p in Path(LIGAND_DIR).rglob("*.sdf")])
    logging.info(f"Found {len(all_files)} SDF files")
    
    if len(all_files) == 0:
        logging.error("No SDF files found")
        return False
    
    # Test 3: File access
    accessible = quick_file_test(all_files[:50])
    logging.info(f"File access test: {len(accessible)}/50 accessible")
    
    # Test 4: Test first few files individually
    logging.info("Testing first 5 files individually...")
    for i, test_file in enumerate(all_files[:5]):
        logging.info(f"Testing file {i+1}/5: {os.path.basename(test_file)}")
        
        test_list = os.path.join(OUTPUT_DIR, f"test_single_{i}.txt")
        with open(test_list, "w") as f:
            f.write(f"{os.path.abspath(test_file)}\n")
        
        # Use the same command structure as chunks
        cmd = ["unidocktools", "mcdock"]
        for flag, value in MCDOCK_FLAGS.items():
            if flag == "--gen_conf":
                cmd.append(flag)
            elif value is not None:
                cmd.extend([flag, value])
        cmd.extend(["--ligand_index", test_list])
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=1800, stdin=subprocess.DEVNULL)
            logging.info(f"  ✓ {os.path.basename(test_file)}: PASSED")
        except subprocess.CalledProcessError as e:
            logging.error(f"  ✗ {os.path.basename(test_file)}: FAILED")
            if "Bad input file" in str(e.stderr):
                logging.error(f"    Bad input file error in {os.path.basename(test_file)}")
            logging.error(f"    STDERR: {e.stderr}")
        except Exception as e:
            logging.error(f"  ✗ {os.path.basename(test_file)}: ERROR - {e}")
        
        # Clean up test file
        if os.path.exists(test_list):
            os.remove(test_list)
    
    logging.info("Diagnostic tests complete")
    return True

def main():
    # Set up signal handlers
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    
    # Handle command line args
    if len(sys.argv) > 1:
        if sys.argv[1] == "--reset":
            reset_progress()
            exit(0)
        elif sys.argv[1] == "--test":
            run_diagnostic_tests()
            exit(0)
        elif sys.argv[1] in ["--help", "-h"]:
            print("Usage: python mcdockv2_b2.py [--reset|--test]")
            exit(0)

    # Validate inputs
    if not os.path.exists(RECEPTOR_FILE):
        logging.error(f"Receptor not found: {RECEPTOR_FILE}")
        exit(1)
    if not os.path.exists(LIGAND_DIR):
        logging.error(f"Ligand dir not found: {LIGAND_DIR}")
        exit(1)
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Collect ligand files
    all_ligand_files = sorted([str(p) for p in Path(LIGAND_DIR).rglob("*.sdf")])
    if not all_ligand_files:
        logging.error("No SDF files found")
        exit(1)

    # Check for completed ligands
    logging.info(f"Found {len(all_ligand_files)} total ligand files")
    completed_ligands = get_completed_ligands()
    
    if completed_ligands:
        logging.info(f"Resume: {len(completed_ligands)} completed")
        ligand_files = filter_completed_ligands(all_ligand_files, completed_ligands)
        logging.info(f"Remaining: {len(ligand_files)}")
        
        if not ligand_files:
            logging.info("All ligands completed")
            exit(0)
    else:
        ligand_files = all_ligand_files

    # Create chunks
    chunks = create_chunks(ligand_files, LIGANDS_PER_CHUNK)
    total_chunks = len(chunks)
    
    logging.info(f"Processing {len(ligand_files)} ligands in {total_chunks} chunks")
    
    # Process chunks
    successful_chunks = 0
    failed_chunks = 0
    
    try:
        for chunk_num, chunk_ligands in enumerate(chunks, 1):
            if shutdown_requested:
                logging.warning("Shutdown requested")
                break
            
            # Check for completed ligands in chunk
            chunk_completed = get_completed_ligands()
            remaining = filter_completed_ligands(chunk_ligands, chunk_completed)
            
            if not remaining:
                logging.info(f"Chunk {chunk_num}: already completed")
                successful_chunks += 1
                continue
            
            # Run chunk
            if run_mcdock_chunk(remaining, chunk_num, total_chunks):
                successful_chunks += 1
                current_completed = get_completed_ligands()
                logging.info(f"Progress: {len(current_completed)}/{len(all_ligand_files)}")
            else:
                failed_chunks += 1
                logging.error(f"Chunk {chunk_num} failed - stopping")
                break
                
    except KeyboardInterrupt:
        logging.error("Interrupted by user")
        current_completed = get_completed_ligands()
        logging.info(f"Progress saved: {len(current_completed)} completed")
        exit(1)
    
    # Final summary
    final_completed = get_completed_ligands()
    total_completed = len(final_completed)
    
    logging.info(f"=== SUMMARY ===")
    logging.info(f"Chunks: {successful_chunks}/{total_chunks} successful")
    logging.info(f"Ligands: {total_completed}/{len(all_ligand_files)} completed")
    
    if failed_chunks > 0:
        logging.error(f"Failed chunks: {failed_chunks}")
        exit(1)
    elif total_completed == len(all_ligand_files):
        logging.info("All ligands completed successfully!")
    else:
        remaining = len(all_ligand_files) - total_completed
        logging.info(f"Remaining: {remaining}")

if __name__ == "__main__":
    main() 
