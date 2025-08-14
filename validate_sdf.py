#!/usr/bin/env python3
"""
Standalone SDF validation script to check for corrupted/malformed SDF files.
Validates all SDF files in a directory using Open Babel and reports issues.
"""
import os
import subprocess
import sys
from pathlib import Path
import logging
import time

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# --- Configuration ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LIGAND_DIR = os.path.join(SCRIPT_DIR, "../data/column_three/ligands_sdf_split/")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "../results/validation/")

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

def validate_all_sdf_files(ligand_dir, output_dir):
    """
    Validate all SDF files in the ligand directory and its subdirectories.
    
    Args:
        ligand_dir (str): Directory containing SDF files to validate
        output_dir (str): Directory to save validation results
        
    Returns:
        tuple: (valid_files, invalid_files, validation_stats)
    """
    if not os.path.exists(ligand_dir):
        logging.error(f"Ligand directory not found: {ligand_dir}")
        return [], [], {}
    
    # Find all SDF files recursively
    all_sdf_files = sorted([str(p) for p in Path(ligand_dir).rglob("*.sdf")])
    
    if not all_sdf_files:
        logging.error(f"No SDF files found in {ligand_dir}")
        return [], [], {}
    
    logging.info(f"Found {len(all_sdf_files)} SDF files to validate")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Validation results
    valid_files = []
    invalid_files = []
    validation_stats = {
        'total_files': len(all_sdf_files),
        'valid_files': 0,
        'invalid_files': 0,
        'start_time': time.time(),
        'validation_errors': {}
    }
    
    # Validate files
    for i, sdf_file in enumerate(all_sdf_files, 1):
        if i % 100 == 0:
            elapsed = time.time() - validation_stats['start_time']
            rate = i / elapsed if elapsed > 0 else 0
            logging.info(f"Validated {i}/{len(all_sdf_files)} files... ({rate:.1f} files/sec)")
        
        if validate_sdf_file(sdf_file):
            valid_files.append(sdf_file)
            validation_stats['valid_files'] += 1
        else:
            invalid_files.append(sdf_file)
            validation_stats['invalid_files'] += 1
            
            # Log the invalid file
            logging.warning(f"Invalid SDF file detected: {os.path.basename(sdf_file)}")
    
    # Calculate final statistics
    validation_stats['elapsed_time'] = time.time() - validation_stats['start_time']
    validation_stats['success_rate'] = (validation_stats['valid_files'] / validation_stats['total_files']) * 100
    
    return valid_files, invalid_files, validation_stats

def save_validation_results(valid_files, invalid_files, validation_stats, output_dir):
    """
    Save validation results to files for later analysis.
    
    Args:
        valid_files (list): List of valid SDF file paths
        invalid_files (list): List of invalid SDF file paths
        validation_stats (dict): Validation statistics
        output_dir (str): Directory to save results
    """
    # Save valid files list
    valid_file_path = os.path.join(output_dir, "valid_sdf_files.txt")
    with open(valid_file_path, 'w') as f:
        for file_path in valid_files:
            f.write(f"{file_path}\n")
    
    # Save invalid files list
    invalid_file_path = os.path.join(output_dir, "invalid_sdf_files.txt")
    with open(invalid_file_path, 'w') as f:
        for file_path in invalid_files:
            f.write(f"{file_path}\n")
    
    # Save validation summary
    summary_path = os.path.join(output_dir, "validation_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("=== SDF VALIDATION SUMMARY ===\n")
        f.write(f"Total files processed: {validation_stats['total_files']}\n")
        f.write(f"Valid files: {validation_stats['valid_files']}\n")
        f.write(f"Invalid files: {validation_stats['invalid_files']}\n")
        f.write(f"Success rate: {validation_stats['success_rate']:.2f}%\n")
        f.write(f"Validation time: {validation_stats['elapsed_time']:.2f} seconds\n")
        f.write(f"Processing rate: {validation_stats['total_files']/validation_stats['elapsed_time']:.1f} files/sec\n")
        
        if invalid_files:
            f.write(f"\n=== INVALID FILES ({len(invalid_files)}) ===\n")
            for invalid_file in invalid_files:
                f.write(f"{os.path.basename(invalid_file)}\n")
    
    logging.info(f"Validation results saved to: {output_dir}")

def main():
    """Main validation function."""
    logging.info("=== SDF VALIDATION SCRIPT ===")
    logging.info(f"Ligand directory: {LIGAND_DIR}")
    logging.info(f"Output directory: {OUTPUT_DIR}")
    
    # Validate all SDF files
    valid_files, invalid_files, validation_stats = validate_all_sdf_files(LIGAND_DIR, OUTPUT_DIR)
    
    if not valid_files and not invalid_files:
        logging.error("No files were processed. Check the ligand directory path.")
        exit(1)
    
    # Save results
    save_validation_results(valid_files, invalid_files, validation_stats, OUTPUT_DIR)
    
    # Print summary
    logging.info(f"\n=== VALIDATION COMPLETE ===")
    logging.info(f"Total files: {validation_stats['total_files']}")
    logging.info(f"Valid files: {validation_stats['valid_files']}")
    logging.info(f"Invalid files: {validation_stats['invalid_files']}")
    logging.info(f"Success rate: {validation_stats['success_rate']:.2f}%")
    logging.info(f"Validation time: {validation_stats['elapsed_time']:.2f} seconds")
    
    if invalid_files:
        logging.warning(f"\nInvalid files found: {len(invalid_files)}")
        logging.warning("Check the validation results for details.")
        logging.warning("Consider re-downloading or fixing invalid files before docking.")
        
        # Show first few invalid files
        logging.warning("First 10 invalid files:")
        for invalid_file in invalid_files[:10]:
            logging.warning(f"  - {os.path.basename(invalid_file)}")
        if len(invalid_files) > 10:
            logging.warning(f"  ... and {len(invalid_files) - 10} more")
    else:
        logging.info("All SDF files are valid! Ready for docking.")
    
    logging.info(f"\nResults saved to: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
