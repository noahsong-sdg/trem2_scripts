#!/usr/bin/env python3
"""
UniDock multi-conformer docking script with comprehensive timing functionality and tranche-aware processing.
This script performs molecular docking using UniDock mcdock and tracks performance metrics across tranches.
"""
import os
import subprocess
from pathlib import Path
import glob

# Import timing utilities
from timing_utils import TimingTracker

# Get the absolute path of the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# --- Configuration ---
RECEPTOR_FILE = os.path.join(SCRIPT_DIR, "../data/receptor/cluster1_receptor.pdbqt")  # PDBQT receptor file for mcdock
LIGAND_DIR = os.path.join(SCRIPT_DIR, "../data/test/ligands_pdbqt_split/")  # Directory containing tranche subdirectories with sdf ligand files
DOCKING_OUTPUT_DIR = os.path.join(SCRIPT_DIR, "../results/mcdock_outputs_test/")

def discover_tranches(ligand_base_dir):
    """
    Discover all tranche directories and count ligands for HTS planning.
    
    Args:
        ligand_base_dir (str): Base directory containing tranche subdirectories
    
    Returns:
        tuple: (tranche_list, total_ligand_count, tranche_stats)
    """
    if not os.path.exists(ligand_base_dir):
        print(f"Error: Ligand directory not found: {ligand_base_dir}")
        return [], 0, {}
    
    tranche_dirs = []
    total_ligands = 0
    tranche_stats = {}
    
    # Find all subdirectories that contain sdf files
    for item in os.listdir(ligand_base_dir):
        item_path = os.path.join(ligand_base_dir, item)
        if os.path.isdir(item_path):
            # Count sdf files in this tranche
            sdf_files = list(Path(item_path).glob("*.sdf"))
            if sdf_files:
                ligand_count = len(sdf_files)
                tranche_dirs.append(item_path)
                tranche_stats[item] = ligand_count
                total_ligands += ligand_count
    
    return tranche_dirs, total_ligands, tranche_stats

def run_tranche_batch_docking(receptor_file, tranche_paths, output_dir, 
                              center_x, center_y, center_z, size_x, size_y, size_z,
                              batch_size=1000, max_tranches_per_batch=5,
                              gen_conf=True, max_num_confs_per_ligand=200, min_rmsd=0.3,
                              scoring_function_rd="vina", exhaustiveness_rd=128, num_modes_rd=3, topn_rd=100,
                              scoring_function_lr="vina", exhaustiveness_lr=512,
                              num_modes_lr=1, topn_lr=1, timer=None):
    """
    Runs UniDock mcdock for multiple tranches with intelligent batching for HTS workflows.

    Args:
        receptor_file (str): Path to the prepared receptor PDB file.
        tranche_paths (list): List of paths to tranche directories containing sdf ligands.
        output_dir (str): Directory to save docking results.
        center_x, center_y, center_z (float): Coordinates of the search space center.
        size_x, size_y, size_z (float): Dimensions of the search space.
        batch_size (int): Maximum number of ligands per batch.
        max_tranches_per_batch (int): Maximum number of tranches to process in one batch.
        gen_conf (bool): Whether to generate conformers for the ligands.
        max_num_confs_per_ligand (int): Maximum number of conformers per ligand.
        min_rmsd (float): Minimum RMSD for conformer generation.
        scoring_function_rd (str): Scoring function for rigid docking.
        exhaustiveness_rd (int): Exhaustiveness for rigid docking.
        num_modes_rd (int): Number of modes for rigid docking.
        topn_rd (int): Top N results for rigid docking.
        scoring_function_lr (str): Scoring function for local refine.
        exhaustiveness_lr (int): Exhaustiveness for local refine.
        num_modes_lr (int): Number of modes for local refine.
        topn_lr (int): Top N results for local refine.
        timer (TimingTracker): Timing tracker instance.
    
    Returns:
        tuple: (successful_dockings, failed_dockings, processed_tranches)
    """
    os.makedirs(output_dir, exist_ok=True)
    
    successful_dockings = 0
    failed_dockings = 0
    processed_tranches = 0
    
    if timer:
        total_ligands = sum(len(list(Path(tp).glob("*.sdf"))) for tp in tranche_paths)
        timer.start_step("Tranche-aware multi-conformer molecular docking", total_ligands)
    
    print(f"🔄 BATCH PROCESSING: {len(tranche_paths)} tranches | Target batch size: {batch_size:,} ligands\n")
    
    # Process tranches in batches
    current_batch_ligands = []
    current_batch_tranches = []
    current_ligand_count = 0
    batch_number = 1
    
    for tranche_path in tranche_paths:
        tranche_name = os.path.basename(tranche_path)
        
        # Get all ligands in this tranche
        ligand_files = list(Path(tranche_path).glob("*.sdf"))
        
        if not ligand_files:
            print(f"⚠️  Warning: No sdf files found in tranche {tranche_name}")
            continue
        
        # Check if adding this tranche exceeds batch limits
        if (len(current_batch_tranches) >= max_tranches_per_batch or 
            current_ligand_count + len(ligand_files) > batch_size):
            
            # Process current batch if it has ligands
            if current_batch_ligands:
                batch_success, batch_failed = _process_ligand_batch(
                    receptor_file, current_batch_ligands, output_dir, batch_number,
                    center_x, center_y, center_z, size_x, size_y, size_z,
                    gen_conf, max_num_confs_per_ligand, min_rmsd,
                    scoring_function_rd, exhaustiveness_rd, num_modes_rd, topn_rd,
                    scoring_function_lr, exhaustiveness_lr, num_modes_lr, topn_lr,
                    timer)
                
                successful_dockings += batch_success
                failed_dockings += batch_failed
                processed_tranches += len(current_batch_tranches)
                
                print(f"🎯 Batch {batch_number} COMPLETE: ✓{batch_success:,} success | ✗{batch_failed:,} failed | 📦{len(current_batch_tranches)} tranches")
                batch_number += 1
            
            # Reset for new batch
            current_batch_ligands = []
            current_batch_tranches = []
            current_ligand_count = 0
        
        # Add tranche to current batch
        current_batch_ligands.extend([str(lf) for lf in ligand_files])
        current_batch_tranches.append(tranche_name)
        current_ligand_count += len(ligand_files)
        
        print(f"  📋 Added tranche {tranche_name}: {len(ligand_files):,} ligands → batch {batch_number}")
    
    # Process final batch
    if current_batch_ligands:
        batch_success, batch_failed = _process_ligand_batch(
            receptor_file, current_batch_ligands, output_dir, batch_number,
            center_x, center_y, center_z, size_x, size_y, size_z,
            gen_conf, max_num_confs_per_ligand, min_rmsd,
            scoring_function_rd, exhaustiveness_rd, num_modes_rd, topn_rd,
            scoring_function_lr, exhaustiveness_lr, num_modes_lr, topn_lr,
            timer)
        
        successful_dockings += batch_success
        failed_dockings += batch_failed
        processed_tranches += len(current_batch_tranches)
        
        print(f"🎯 Batch {batch_number} COMPLETE: ✓{batch_success:,} success | ✗{batch_failed:,} failed | 📦{len(current_batch_tranches)} tranches")
    
    if timer:
        timer.end_step()
    
    return successful_dockings, failed_dockings, processed_tranches

def _process_ligand_batch(receptor_file, ligand_files, output_dir, batch_number,
                         center_x, center_y, center_z, size_x, size_y, size_z,
                         gen_conf, max_num_confs_per_ligand, min_rmsd,
                         scoring_function_rd, exhaustiveness_rd, num_modes_rd, topn_rd,
                         scoring_function_lr, exhaustiveness_lr, num_modes_lr, topn_lr,
                         timer):
    """Process a batch of ligands through mcdock."""
    
    # Validate ligand file format
    sample_file = ligand_files[0] if ligand_files else ""
    if sample_file.endswith('.pdbqt'):
        print(f"⚠️  WARNING: mcdock expects SDF format, but ligands are in pdbqt format")
        print(f"   This may cause compatibility issues. Consider converting to SDF format.")
    
    # Create batch-specific output directory
    batch_output_dir = os.path.join(output_dir, f"batch_{batch_number}")
    os.makedirs(batch_output_dir, exist_ok=True)
    
    # Write ligand list to a file to avoid "Argument list too long" error
    ligand_list_file = os.path.join(batch_output_dir, f"ligand_list_batch_{batch_number}.txt")
    with open(ligand_list_file, 'w') as f:
        for ligand_file in ligand_files:
            f.write(f"{os.path.abspath(ligand_file)}\n")
    
    print(f"📝 Created ligand list file: {ligand_list_file} ({len(ligand_files):,} ligands)")
    
    # Build mcdock command using ligand file input
    command = [
        "unidocktools", "mcdock",
        "--receptor", os.path.abspath(receptor_file),
        "--ligand_index", ligand_list_file,  # Use file input instead of comma-separated list
        "--center_x", str(center_x),
        "--center_y", str(center_y), 
        "--center_z", str(center_z),
        "--size_x", str(size_x),
        "--size_y", str(size_y),
        "--size_z", str(size_z),
        "--workdir", os.path.join(batch_output_dir, "MultiConfDock"),
        "--savedir", os.path.join(batch_output_dir, "MultiConfDock-Result"),
        "--batch_size", str(min(len(ligand_files), 1200)),  # Adjust batch size for available ligands
        "--scoring_function_rigid_docking", scoring_function_rd,
        "--exhaustiveness_rigid_docking", str(exhaustiveness_rd),
        "--num_modes_rigid_docking", str(num_modes_rd),
        "--topn_rigid_docking", str(topn_rd),
        "--scoring_function_local_refine", scoring_function_lr,
        "--exhaustiveness_local_refine", str(exhaustiveness_lr),
        "--num_modes_local_refine", str(num_modes_lr),
        "--topn_local_refine", str(topn_lr),
        "--min_rmsd", str(min_rmsd),
        "--max_num_confs_per_ligand", str(max_num_confs_per_ligand)
    ]
    
    # Add conformer generation flag if requested
    if gen_conf:
        command.append("--gen_conf")
    
    print(f"\n{'='*60}")
    print(f"🧪 BATCH {batch_number}: Processing {len(ligand_files):,} ligands")
    print(f"📁 Output: batch_{batch_number}/")
    print(f"⚙️  UniDock mcdock starting...")
    print(f"{'='*60}")
    
    successful_dockings = 0
    failed_dockings = 0
    
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        
        # Check for successful outputs in the result directory
        result_dir = os.path.join(batch_output_dir, "MultiConfDock-Result")
        if os.path.exists(result_dir):
            # Count result files
            result_files = []
            for root, dirs, files in os.walk(result_dir):
                for file in files:
                    if file.endswith(('.sdf', '.sdf', '.csv')):
                        result_files.append(file)
            
            if result_files:
                successful_dockings = len(ligand_files)  # Assume all succeeded if results exist
                print(f"🎉 Batch {batch_number} SUCCESS: {len(result_files):,} result files generated")
            else:
                failed_dockings = len(ligand_files)
                print(f"❌ Batch {batch_number} FAILED: No result files found")
        else:
            failed_dockings = len(ligand_files)
            print(f"❌ Batch {batch_number} FAILED: Result directory missing")
        
        if timer:
            timer.update_progress(len(ligand_files))
        
    except subprocess.CalledProcessError as e:
        print(f"✗ Error during batch {batch_number} UniDock mcdock execution:")
        print(f"Return code: {e.returncode}")
        if e.stdout:
            print(f"Standard output: {e.stdout[:500]}...")
        if e.stderr:
            print(f"Error output: {e.stderr[:500]}...")  # Truncate long error messages
        failed_dockings = len(ligand_files)
        
    except FileNotFoundError:
        print(f"✗ Error: unidocktools not found. Please ensure it's in your PATH.")
        failed_dockings = len(ligand_files)
    except Exception as e:
        print(f"✗ Unexpected error during batch {batch_number}: {e}")
        failed_dockings = len(ligand_files)
    
    return successful_dockings, failed_dockings

# --- Run UniDock mcdock ---
def run_unidock_mcdock(receptor_file, ligand_input, output_dir, 
                       center_x, center_y, center_z, size_x, size_y, size_z,
                       gen_conf=True, max_num_confs_per_ligand=200, min_rmsd=0.3,
                       batch_size=1200, scoring_function_rd="vina", 
                       exhaustiveness_rd=128, num_modes_rd=3, topn_rd=100,
                       scoring_function_lr="vina", exhaustiveness_lr=512,
                       num_modes_lr=1, topn_lr=1, timer=None):
    """
    Runs UniDock mcdock for multi-conformer docking with automatic tranche detection and batch processing.

    Args:
        receptor_file (str): Path to the prepared receptor PDBQT file (mcdock requires PDBQT format).
        ligand_input (str): Path to the directory containing tranche subdirectories with sdf ligand files 
                           OR path to a single sdf ligand file.
        output_dir (str): Directory to save docking results.
        center_x, center_y, center_z (float): Coordinates of the search space center.
        size_x, size_y, size_z (float): Dimensions of the search space (Angstroms).
        gen_conf (bool): Whether to generate conformers for the ligands.
        max_num_confs_per_ligand (int): Maximum number of conformers to generate per ligand.
        min_rmsd (float): Minimum RMSD for conformer generation.
        batch_size (int): Batch size for docking.
        scoring_function_rd (str): Scoring function for rigid docking.
        exhaustiveness_rd (int): Exhaustiveness for rigid docking.
        num_modes_rd (int): Number of modes for rigid docking.
        topn_rd (int): Top N results for rigid docking.
        scoring_function_lr (str): Scoring function for local refine.
        exhaustiveness_lr (int): Exhaustiveness for local refine.
        num_modes_lr (int): Number of modes for local refine.
        topn_lr (int): Top N results for local refine.
        timer (TimingTracker): Timing tracker instance.
    
    Returns:
        tuple: (successful_dockings, failed_dockings, processed_tranches)
    """
    os.makedirs(output_dir, exist_ok=True)

    # Check if this is a tranche-based directory structure or single file/flat directory
    if os.path.isfile(ligand_input) and ligand_input.endswith('.sdf'):
        # Single ligand file - use original logic
        print("Processing single ligand file...")
        ligand_files = [ligand_input]
        successful, failed = _process_ligand_batch(
            receptor_file, ligand_files, output_dir, 1,
            center_x, center_y, center_z, size_x, size_y, size_z,
            gen_conf, max_num_confs_per_ligand, min_rmsd,
            scoring_function_rd, exhaustiveness_rd, num_modes_rd, topn_rd,
            scoring_function_lr, exhaustiveness_lr, num_modes_lr, topn_lr,
            timer)
        return successful, failed, 1
    
    elif os.path.isdir(ligand_input):
        # Check for tranche-based structure
        tranche_dirs, total_ligands, tranche_stats = discover_tranches(ligand_input)
        
        if tranche_dirs:
            # Tranche-based structure detected
            print(f"\n🔍 TRANCHE DISCOVERY")
            print(f"📊 Found {len(tranche_dirs)} tranches with {total_ligands:,} total ligands")
            print(f"📋 Preparing batch processing...\n")
            
            return run_tranche_batch_docking(
                receptor_file, tranche_dirs, output_dir,
                center_x, center_y, center_z, size_x, size_y, size_z,
                batch_size=batch_size, max_tranches_per_batch=3,  # Reduced to prevent arg list issues
                gen_conf=gen_conf, max_num_confs_per_ligand=max_num_confs_per_ligand,
                min_rmsd=min_rmsd, scoring_function_rd=scoring_function_rd,
                exhaustiveness_rd=exhaustiveness_rd, num_modes_rd=num_modes_rd,
                topn_rd=topn_rd, scoring_function_lr=scoring_function_lr,
                exhaustiveness_lr=exhaustiveness_lr, num_modes_lr=num_modes_lr,
                topn_lr=topn_lr, timer=timer)
        else:
            # Flat directory structure - collect all sdf files
            ligand_files = list(Path(ligand_input).glob("*.sdf"))
            if ligand_files:
                print(f"Processing flat directory with {len(ligand_files)} ligands...")
                ligand_file_paths = [str(lf) for lf in ligand_files]
                successful, failed = _process_ligand_batch(
                    receptor_file, ligand_file_paths, output_dir, 1,
                    center_x, center_y, center_z, size_x, size_y, size_z,
                    gen_conf, max_num_confs_per_ligand, min_rmsd,
                    scoring_function_rd, exhaustiveness_rd, num_modes_rd, topn_rd,
                    scoring_function_lr, exhaustiveness_lr, num_modes_lr, topn_lr,
                    timer)
                return successful, failed, 1
    
    print(f"No valid sdf ligand files found in {ligand_input}")
    return 0, 0, 0

if __name__ == "__main__":
    # Initialize timing tracker
    timer = TimingTracker("03_mcdock_tranche_aware")
    
    try:
        # --- !!! IMPORTANT: Define Receptor Binding Site !!! ---
        # These are placeholder values. Replace with actual coordinates for TREM2.
        # You can get these from literature, by analyzing the receptor structure in a
        # molecular viewer (like PyMOL, ChimeraX) with a known ligand, or using pocket detection tools.
        CENTER_X, CENTER_Y, CENTER_Z = 42.328, 28.604, 21.648 
        SIZE_X, SIZE_Y, SIZE_Z = 22.5, 22.5, 22.5 # Default mcdock box size

        print("=== UniDock Multi-Conformer HTS Workflow with Tranche-Aware Processing ===")

        # 1. Check for receptor and ligands
        print("\n--- Step 1: Validating Input Files ---")
        timer.start_step("Validate input files")
        
        if not os.path.exists(RECEPTOR_FILE):
            print(f"\nError: Receptor file not found at {RECEPTOR_FILE}")
            print("Please prepare your TREM2 receptor in PDBQT format and place it there.")
            print("Note: mcdock requires PDBQT format, not SDF.")
            timer.finish()
            exit(1)

        if not os.path.exists(LIGAND_DIR):
            print(f"\nError: Ligand directory not found at {LIGAND_DIR}")
            print("Please ensure the tranche-organized sdf ligand files exist.")
            timer.finish()
            exit(1)
        
        # Discover tranches and validate structure
        tranche_dirs, total_ligands, tranche_stats = discover_tranches(LIGAND_DIR)
        
        if not tranche_dirs:
            print(f"\nError: No valid tranche directories with sdf files found in {LIGAND_DIR}")
            print("Expected directory structure: ligands_sdf_split/TRANCHE_NAME/*.sdf")
            timer.finish()
            exit(1)
        
        print(f"\n✓ Discovered {len(tranche_dirs)} tranches with {total_ligands} total ligands")
        print("Tranche overview:")
        for tranche_name, count in sorted(tranche_stats.items()):
            print(f"  {tranche_name}: {count:,} ligands")
        
        timer.end_step()
        
        # 2. Run UniDock mcdock with tranche-aware processing
        print("\n--- Step 2: Running Tranche-Aware UniDock Multi-Conformer Docking ---")
        print("Processing ligands by tranche with intelligent batching for optimal HTS performance.")
        
        successful_dockings, failed_dockings, processed_tranches = run_unidock_mcdock(
            RECEPTOR_FILE, LIGAND_DIR, DOCKING_OUTPUT_DIR, 
            CENTER_X, CENTER_Y, CENTER_Z, SIZE_X, SIZE_Y, SIZE_Z,
            gen_conf=True,
            max_num_confs_per_ligand=50,  # Reduced for first screening
            min_rmsd=0.3,  # Slightly higher for faster generation
            batch_size=1200,
            scoring_function_rd="vina",  # Rigid docking scoring
            exhaustiveness_rd=32,  # Reduced for first screening speed
            num_modes_rd=3,  # Rigid docking modes
            topn_rd=20,  # Fewer candidates for refinement
            scoring_function_lr="vina",  # Local refine scoring
            exhaustiveness_lr=64,  # Much lower for screening
            num_modes_lr=1,  # Local refine modes
            topn_lr=1,  # Final top N results
            timer=timer,
            
        )
        
        # Set final ligand count for performance metrics
        timer.set_final_ligand_count(successful_dockings)
        
        # Generate final timing report
        report = timer.finish()
        
        print(f"\n=== TRANCHE-AWARE MULTI-CONFORMER DOCKING SUMMARY ===")
        print(f"✓ Successful dockings: {successful_dockings:,}")
        print(f"✗ Failed dockings: {failed_dockings:,}")
        print(f"📁 Tranches processed: {processed_tranches}")
        print(f"📁 Docking outputs saved to: {DOCKING_OUTPUT_DIR}")
        
        # Performance summary
        if "performance_metrics" in report:
            metrics = report["performance_metrics"]
            print(f"\n📊 Tranche-Aware HTS Performance Metrics:")
            print(f"   Docking rate: {metrics['ligands_per_minute']:.1f} ligands/minute")
            print(f"   Average time per ligand: {metrics['average_seconds_per_ligand']:.3f} seconds")
            print(f"   Tranches per hour: {(processed_tranches / (report['total_time_seconds'] / 3600)):.1f}")
            print(f"\n🚀 HTS Scale Estimates:")
            print(f"   1M ligands would take: {metrics['estimated_time_for_1M_ligands']}")
            print(f"   10M ligands would take: {metrics['estimated_time_for_10M_ligands']}")
            
            # Tranche-specific estimates
            if processed_tranches > 0:
                avg_ligands_per_tranche = successful_dockings / processed_tranches
                print(f"   Average ligands per tranche: {avg_ligands_per_tranche:.0f}")
    
    except Exception as e:
        print(f"Error during tranche-aware multi-conformer docking workflow: {e}")
        timer.finish()
        exit(1) 
