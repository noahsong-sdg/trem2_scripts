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

def run_mcdock_chunk_with_retry(chunk_ligands, chunk_num, total_chunks):
    """Run UniDock on a chunk of ligands, skipping failed individual ligands."""
    logging.info(f"Chunk {chunk_num}/{total_chunks}: {len(chunk_ligands)} ligands")
    
    # Quick file test
    accessible_files = quick_file_test(chunk_ligands)
    if len(accessible_files) != len(chunk_ligands):
        logging.warning(f"Chunk {chunk_num}: {len(chunk_ligands) - len(accessible_files)} inaccessible files")
        chunk_ligands = accessible_files  # Use only accessible files
    
    if not chunk_ligands:
        logging.error(f"Chunk {chunk_num}: No accessible files")
        return False
    
    # Track failed ligands
    failed_ligands = []
    successful_ligands = []
    remaining_ligands = chunk_ligands.copy()
    
    # Try processing the chunk multiple times, removing failed ligands each time
    max_attempts = 3
    for attempt in range(max_attempts):
        if not remaining_ligands:
            break
            
        logging.info(f"Chunk {chunk_num} attempt {attempt + 1}: {len(remaining_ligands)} ligands")
        
        # Create ligand list file for remaining ligands
        chunk_list = os.path.join(OUTPUT_DIR, f"chunk_{chunk_num}_attempt_{attempt + 1}_ligand_list.txt")
        with open(chunk_list, "w") as f:
            for ligand in remaining_ligands:
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
            logging.info(f"Chunk {chunk_num} attempt {attempt + 1} completed successfully")
            successful_ligands.extend(remaining_ligands)
            os.remove(chunk_list)
            break
            
        except subprocess.CalledProcessError as e:
            logging.warning(f"Chunk {chunk_num} attempt {attempt + 1} failed")
            
            # Extract failed ligands from error message
            failed_in_attempt = extract_failed_ligands_from_error(e.stderr, remaining_ligands)
            
            if failed_in_attempt:
                logging.warning(f"Identified {len(failed_in_attempt)} failed ligands in attempt {attempt + 1}")
                failed_ligands.extend(failed_in_attempt)
                
                # Remove failed ligands from remaining list
                remaining_ligands = [ligand for ligand in remaining_ligands if ligand not in failed_in_attempt]
                
                if remaining_ligands:
                    logging.info(f"Retrying with {len(remaining_ligands)} remaining ligands")
                else:
                    logging.error(f"No ligands remaining after removing failed ones")
                    break
            else:
                # Couldn't identify specific failed ligands, try a different approach
                logging.warning(f"Chunk {chunk_num} attempt {attempt + 1}: Couldn't identify specific failed ligands")
                logging.warning("This might be a different type of error. Logging details:")
                if e.stderr:
                    logging.warning(f"STDERR: {e.stderr}")
                if e.stdout:
                    logging.warning(f"STDOUT: {e.stdout}")
                
                # If this is the last attempt, mark all remaining ligands as failed
                if attempt == max_attempts - 1:
                    logging.error(f"Marking all {len(remaining_ligands)} remaining ligands as failed")
                    failed_ligands.extend(remaining_ligands)
                    remaining_ligands = []
                break
                
        except subprocess.TimeoutExpired:
            logging.error(f"Chunk {chunk_num} attempt {attempt + 1} timed out")
            break
        finally:
            # Clean up attempt file
            if os.path.exists(chunk_list):
                os.remove(chunk_list)
    
    # Log results
    logging.info(f"Chunk {chunk_num} results:")
    logging.info(f"  Successful: {len(successful_ligands)} ligands")
    logging.info(f"  Failed: {len(failed_ligands)} ligands")
    
    # Save failed ligands to file
    if failed_ligands:
        failed_file = os.path.join(OUTPUT_DIR, f"chunk_{chunk_num}_failed_ligands.txt")
        with open(failed_file, "w") as f:
            for ligand in failed_ligands:
                f.write(f"{ligand}\n")
        logging.info(f"Failed ligands saved to: {failed_file}")
    
    # Return True if at least some ligands were processed successfully
    success = len(successful_ligands) > 0
    if not success:
        logging.error(f"Chunk {chunk_num}: All ligands failed processing")
        logging.error(f"  Total ligands in chunk: {len(chunk_ligands)}")
        logging.error(f"  Failed ligands: {len(failed_ligands)}")
        logging.error(f"  Successful ligands: {len(successful_ligands)}")
    return success

def extract_failed_ligands_from_error(stderr, ligand_list):
    """Extract list of failed ligands from UniDock error message."""
    if not stderr:
        return []
    
    failed_ligands = []
    error_text = str(stderr)
    
    # Look for "Bad input file" errors
    if "Bad input file" in error_text:
        # Extract filenames from error messages
        error_lines = error_text.split('\n')
        for line in error_lines:
            if "Bad input file" in line and ".sdf" in line:
                # Try to extract the filename
                if "obabel_" in line:
                    parts = line.split('/')
                    for part in parts:
                        if part.endswith('.sdf'):
                            # Find the corresponding ligand file path
                            for ligand in ligand_list:
                                if part in ligand:
                                    failed_ligands.append(ligand)
                                    logging.warning(f"Identified failed ligand: {os.path.basename(ligand)}")
                                    break
                            break
    
    return failed_ligands

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
            
            # Run chunk with retry logic
            chunk_success = run_mcdock_chunk_with_retry(remaining, chunk_num, total_chunks)
            
            if chunk_success:
                successful_chunks += 1
                current_completed = get_completed_ligands()
                logging.info(f"Progress: {len(current_completed)}/{len(all_ligand_files)}")
            else:
                failed_chunks += 1
                logging.warning(f"Chunk {chunk_num} failed - continuing with remaining chunks")
                # Don't break, continue to next chunk
                
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
        logging.warning(f"Failed chunks: {failed_chunks}")
        logging.info("Check failed ligand files in output directory for details")
        logging.info("You can re-run the script to retry failed chunks")
    elif total_completed == len(all_ligand_files):
        logging.info("All ligands completed successfully!")
    else:
        remaining = len(all_ligand_files) - total_completed
        logging.info(f"Remaining: {remaining}")
        logging.info("Run the script again to resume from where it left off")
    
    # Don't exit with error code if some chunks failed but we processed successfully
    if total_completed > 0:
        logging.info("Script completed with partial success")
    else:
        logging.error("No ligands were processed successfully")
        exit(1)

if __name__ == "__main__":
    main() 
