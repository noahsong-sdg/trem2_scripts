#!/usr/bin/env python3
"""
SDF Data download script with comprehensive timing functionality and parallel downloads.
Downloads ligand data in SDF format for mcdock compatibility and tracks performance metrics.
conda create -n unidock_env unidock -c conda-forge
mamba create -n unidock_env python=3.10 requests unidock -c conda-forge'




ctrl f column_three
repalce
"""
import requests
import os
import gzip
import shutil
import time
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

# Get the absolute path of the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Thread-safe progress tracking
download_lock = threading.Lock()
progress_counter = {'completed': 0, 'failed': 0}

def download_zinc_subset(url, output_dir, filename=None):
    """
    Downloads a subset of the ZINC database in SDF format.

    Args:
        url (str): URL to the ZINC subset (e.g., a direct download link to a SDF.gz file).
        output_dir (str): Directory to save the downloaded file.
        filename (str): Name of the file to save. If None, extracts from URL.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Extract filename from URL if not provided
    if filename is None:
        filename = url.split('/')[-1]

    filepath = os.path.join(output_dir, filename)

    # Existence check: skip if file exists and is non-empty
    if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
        logging.info(f"✓ Skipping already downloaded: {filename}")
        return filepath

    try:
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        # Get file size if available
        total_size = int(response.headers.get('content-length', 0))
        downloaded_size = 0
        
        with open(filepath, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                downloaded_size += len(chunk)
        
        file_size = os.path.getsize(filepath)
        
        # Thread-safe progress update
        with download_lock:
            progress_counter['completed'] += 1
            logging.info(f"✓ Downloaded ({progress_counter['completed']}) {filename} ({file_size:,} bytes)")
        
        return filepath
        
    except requests.exceptions.RequestException as e:
        with download_lock:
            progress_counter['failed'] += 1
            logging.error(f"✗ Failed ({progress_counter['failed']}) {filename}: {e}")
        return None
    except Exception as e:
        with download_lock:
            progress_counter['failed'] += 1
            logging.error(f"✗ Error ({progress_counter['failed']}) {filename}: {e}")
        return None

def download_single_file(args):
    """Helper function for parallel downloads"""
    url, output_dir, filename = args
    return download_zinc_subset(url, output_dir, filename)

def download_all_from_uri_file(uri_file_path, base_output_dir, max_workers=4, max_retries=5):
    """
    Reads a list of URLs from a .uri file and downloads SDF data from each using parallel processing.
    Retries failed downloads up to max_retries times.

    Args:
        uri_file_path (str): Path to the .uri file containing URLs.
        base_output_dir (str): The base directory to save downloaded files.
        max_workers (int): Number of parallel download threads.
        max_retries (int): Maximum number of retry attempts for failed downloads.
    """
    if not os.path.exists(uri_file_path):
        logging.error(f"Error: URI file not found at {uri_file_path}")
        return 0, 0

    # Read all URLs
    urls = []
    with open(uri_file_path, 'r') as f:
        for line in f:
            url = line.strip()
            if url and not url.startswith("#"):
                urls.append(url)
    
    if not urls:
        logging.error("No valid URLs found in URI file")
        return 0, 0
    
    logging.info(f"Starting parallel download of {len(urls)} SDF files using {max_workers} workers...")
    logging.info(f"Will retry failed downloads up to {max_retries} times")
    
    # Reset progress counters
    with download_lock:
        progress_counter['completed'] = 0
        progress_counter['failed'] = 0
    
    # Prepare initial download arguments
    download_args = []
    for i, url in enumerate(urls, 1):
        filename = url.split('/')[-1]
        if not filename:
            filename = f"downloaded_ligand_{i}.sdf.gz"
        download_args.append((url, base_output_dir, filename))
    
    # Track failed downloads for retries
    failed_downloads = []
    total_downloaded = 0
    total_failed = 0
    
    # Execute downloads with retries
    for attempt in range(1, max_retries + 1):
        if attempt == 1:
            logging.info(f"\n=== ATTEMPT {attempt}: Initial download ===")
            current_args = download_args
        else:
            if not failed_downloads:
                logging.info(f"\n=== All downloads successful after {attempt-1} attempts ===")
                break
            logging.info(f"\n=== ATTEMPT {attempt}: Retrying {len(failed_downloads)} failed downloads ===")
            current_args = failed_downloads
            failed_downloads = []  # Reset for next attempt
        
        # Reset progress counters for this attempt
        with download_lock:
            progress_counter['completed'] = 0
            progress_counter['failed'] = 0
        
        # Execute parallel downloads for current attempt
        attempt_downloaded = 0
        attempt_failed = 0
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all download tasks for this attempt
            future_to_args = {executor.submit(download_single_file, args): args 
                            for args in current_args}
            
            # Process completed downloads
            for future in as_completed(future_to_args):
                args = future_to_args[future]
                url = args[0]
                try:
                    result = future.result()
                    if result:
                        attempt_downloaded += 1
                    else:
                        attempt_failed += 1
                        failed_downloads.append(args)  # Add to retry list
                except Exception as e:
                    attempt_failed += 1
                    failed_downloads.append(args)  # Add to retry list
                    logging.error(f"✗ Exception downloading {url}: {e}")
        
        total_downloaded += attempt_downloaded
        total_failed = len(failed_downloads)
        
        logging.info(f"Attempt {attempt} results: {attempt_downloaded} successful, {attempt_failed} failed")
        
        if attempt < max_retries and failed_downloads:
            logging.info(f"Retrying {len(failed_downloads)} failed downloads in next attempt...")
            logging.info(f"Waiting 60 minutes before retry attempt {attempt + 1}...")
            time.sleep(3600)  # 60 minutes = 600 seconds
        elif failed_downloads:
            logging.warning(f"Reached maximum retries ({max_retries}). {len(failed_downloads)} downloads still failed.")
    
    logging.info(f"\n=== FINAL DOWNLOAD SUMMARY ===")
    logging.info(f"✓ Total successful downloads: {total_downloaded}")
    logging.info(f"✗ Total failed downloads: {total_failed}")
    logging.info(f"📊 Success rate: {total_downloaded/(total_downloaded+total_failed)*100:.1f}%")
    
    return total_downloaded, total_failed

def extract_sdf_files(raw_dir, output_dir, max_workers=4):
    """
    Extract .sdf.gz files to individual .sdf files ready for docking with parallel processing.
    
    Args:
        raw_dir (str): Directory containing downloaded .sdf.gz files
        output_dir (str): Directory to save extracted .sdf files
        max_workers (int): Number of parallel extraction threads
    
    Returns:
        tuple: (successful_extractions, failed_extractions)
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Find all .sdf.gz files
    gz_files = list(Path(raw_dir).glob("*.sdf.gz"))
    
    if not gz_files:
        logging.warning(f"No .sdf.gz files found in {raw_dir}")
        return 0, 0
    
    logging.info(f"Extracting {len(gz_files)} SDF files using {max_workers} workers...")
    
    def extract_single_file(gz_file):
        """Helper function for parallel extraction"""
        try:
            # Create output filename (remove .gz extension)
            output_file = Path(output_dir) / gz_file.stem
            
            # Extract the gzipped file
            with gzip.open(gz_file, 'rb') as f_in:
                with open(output_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            return gz_file.name, True, None
            
        except Exception as e:
            return gz_file.name, False, str(e)
    
    successful = 0
    failed = 0
    
    # Execute parallel extractions
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all extraction tasks
        future_to_file = {executor.submit(extract_single_file, gz_file): gz_file 
                         for gz_file in gz_files}
        
        # Process completed extractions
        for future in as_completed(future_to_file):
            gz_file = future_to_file[future]
            try:
                filename, success, error = future.result()
                if success:
                    successful += 1
                    logging.info(f"✓ Extracted ({successful}): {filename}")
                else:
                    failed += 1
                    logging.error(f"✗ Failed ({failed}): {filename} - {error}")
            except Exception as e:
                failed += 1
                logging.error(f"✗ Exception extracting {gz_file.name}: {e}")
    
    logging.info(f"Extraction complete! Successful: {successful}, Failed: {failed}")
    return successful, failed

def get_tranche_name_from_filename(filename):
    """Extract tranche name from ZINC filename for organization."""
    # ZINC files are typically named like: BFEDMM.xaa.sdf.gz
    # Extract the tranche part (e.g., BFEDMM.xaa)
    if '.sdf' in filename:
        base_name = filename.replace('.sdf.gz', '').replace('.sdf', '')
        if '.' in base_name:
            return base_name
    return "unknown_tranche"

def _save_molecule(molecule_lines, molecule_name, output_dir, tranche_name, molecule_index):
    """Save a single molecule to its own SDF file in appropriate tranche directory."""
    if not molecule_lines:
        return
    
    # Test if Open Babel can read this molecule
    try:
        import subprocess
        mol_text = '\n'.join(molecule_lines)
        result = subprocess.run(['obabel', '-:', '-osdf'], 
                              input=mol_text, text=True, capture_output=True, timeout=10)
        if result.returncode != 0:
            logging.warning(f"Skipping invalid molecule {molecule_index} in {tranche_name}")
            return
    except subprocess.TimeoutExpired:
        logging.warning(f"Skipping molecule {molecule_index} in {tranche_name} (timeout)")
        return
    except FileNotFoundError:
        # Open Babel not available, continue without validation
        pass
    except Exception as e:
        logging.error(f"Error validating molecule {molecule_index}: {e}")
        return
    
    # Create tranche-specific directory
    tranche_dir = os.path.join(output_dir, tranche_name)
    os.makedirs(tranche_dir, exist_ok=True)
    
    # Generate filename
    if molecule_name:
        # Use molecule name if available (first line in SDF)
        filename = f"{molecule_name}.sdf"
        # Clean filename - remove invalid characters
        filename = "".join(c for c in filename if c.isalnum() or c in "._-")
    else:
        # Extract ZINC ID from molecule content if available
        zinc_id = None
        for line in molecule_lines[:5]:  # Check first few lines for ZINC ID
            if 'ZINC' in line.strip().upper():
                # Extract ZINC ID from the line
                parts = line.strip().split()
                for part in parts:
                    if part.upper().startswith('ZINC'):
                        zinc_id = part
                        break
                if zinc_id:
                    break
        
        if zinc_id:
            filename = f"{zinc_id}.sdf"
        else:
            filename = f"molecule_{molecule_index:06d}.sdf"
    
    output_path = os.path.join(tranche_dir, filename)
    
    with open(output_path, 'w') as f:
        for line in molecule_lines:
            # Filter out unwanted header lines
            """line_stripped = line.strip()
            if (line_stripped.startswith('OpenBabel') or 
                line_stripped.startswith('ZINC') or
                line_stripped.endswith('.sdf') or
                line_stripped.endswith('.pdbqt') or
                (len(line_stripped) < 5 and not line_stripped.startswith('$$$$'))):
                continue"""
            f.write(line + '\n')

def split_sdf_files(input_dir, output_dir, max_workers=4):
    """
    Split multi-molecule SDF files into individual single-molecule files organized by tranche.
    This is required for UniDock mcdock which expects one molecule per file.
    
    Args:
        input_dir (str): Directory containing multi-molecule SDF files
        output_dir (str): Directory to save individual SDF files (organized by tranche)
        max_workers (int): Number of parallel splitting threads
    
    Returns:
        tuple: (total_molecules_split, failed_files, tranche_count)
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Find all .sdf files
    sdf_files = list(Path(input_dir).glob("*.sdf"))
    
    if not sdf_files:
        logging.warning(f"No .sdf files found in {input_dir}")
        return 0, 0, 0
    
    logging.info(f"Splitting {len(sdf_files)} SDF files into tranche-organized individual molecules using {max_workers} workers...")
    
    def split_single_sdf(sdf_file):
        """Helper function for parallel SDF splitting"""
        try:
            molecule_count = 0
            current_molecule = []
            current_name = None
            tranche_name = get_tranche_name_from_filename(sdf_file.name)
            
            with open(sdf_file, 'r') as f:
                for line in f:
                    line = line.rstrip('\n\r')
                    
                    if line.strip() == '$$$$':
                        # End of molecule in SDF format
                        if current_molecule:
                            current_molecule.append(line)
                            _save_molecule(current_molecule, current_name, output_dir, 
                                         tranche_name, molecule_count)
                            molecule_count += 1
                            current_molecule = []
                            current_name = None
                    else:
                        current_molecule.append(line)
                        
                        # Extract molecule name if it's the first line (title line in SDF)
                        if len(current_molecule) == 1 and line.strip():
                            current_name = line.strip()
            
            # Handle case where file doesn't end with $$$$
            if current_molecule:
                current_molecule.append('$$$$')
                _save_molecule(current_molecule, current_name, output_dir, 
                             tranche_name, molecule_count)
                molecule_count += 1
            
            return sdf_file.name, molecule_count, tranche_name, None
            
        except Exception as e:
            return sdf_file.name, 0, "unknown", str(e)
    
    total_molecules = 0
    failed_files = 0
    tranches_created = set()
    
    # Execute parallel splitting
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all splitting tasks
        future_to_file = {executor.submit(split_single_sdf, sdf_file): sdf_file 
                         for sdf_file in sdf_files}
        
        # Process completed splits
        for future in as_completed(future_to_file):
            sdf_file = future_to_file[future]
            try:
                filename, molecule_count, tranche_name, error = future.result()
                if error:
                    failed_files += 1
                    logging.error(f"✗ Failed to split {filename}: {error}")
                else:
                    total_molecules += molecule_count
                    tranches_created.add(tranche_name)
                    logging.info(f"✓ Split {filename}: {molecule_count} molecules → tranche {tranche_name}")
            except Exception as e:
                failed_files += 1
                logging.error(f"✗ Exception splitting {sdf_file.name}: {e}")
    
    logging.info(f"Splitting complete! Total molecules: {total_molecules}, Failed files: {failed_files}, Tranches: {len(tranches_created)}")
    return total_molecules, failed_files, len(tranches_created)

if __name__ == "__main__":
    
    try:
        RAW_LIGANDS_DIR = os.path.join(SCRIPT_DIR, "../data/column_three/ligands_raw")
        URI_FILE = os.path.join(SCRIPT_DIR, "../data/column_three.uri") # Using the SDF.gz file URLs
        
        # Configuration for parallel processing
        DOWNLOAD_WORKERS = 8  # Number of parallel download threads
        EXTRACTION_WORKERS = 4  # Number of parallel extraction threads
        
        logging.info(f"=== PARALLEL SDF DATA DOWNLOAD SCRIPT ===")
        logging.info(f"Download workers: {DOWNLOAD_WORKERS}")
        logging.info(f"Extraction workers: {EXTRACTION_WORKERS}")
        
        # Check if URI file exists
        if not os.path.exists(URI_FILE):
            logging.error(f"Error: URI file not found at {URI_FILE}")
            logging.error("Please ensure the URI file exists with URLs to SDF downloads.")
            exit(1)
        
        # Download all files with parallel processing and retries
        successful_downloads, failed_downloads = download_all_from_uri_file(
            URI_FILE, RAW_LIGANDS_DIR, max_workers=DOWNLOAD_WORKERS, max_retries=5)
                
        logging.info(f"\n=== DOWNLOAD SUMMARY ===")
        logging.info(f"✓ Successful downloads: {successful_downloads}")
        logging.info(f"✗ Failed downloads: {failed_downloads}")
        logging.info(f"📁 Files saved to: {RAW_LIGANDS_DIR}")
        
        # Extract SDF files with parallel processing
        sdf_dir = os.path.join(SCRIPT_DIR, "../data/column_three/ligands_sdf")
        successful_extractions, failed_extractions = extract_sdf_files(
            RAW_LIGANDS_DIR, sdf_dir, max_workers=EXTRACTION_WORKERS)
        logging.info(f"\n=== EXTRACTION SUMMARY ===")
        logging.info(f"✓ Successful extractions: {successful_extractions}")
        logging.info(f"✗ Failed extractions: {failed_extractions}")
        logging.info(f"📁 SDF files ready for splitting: {sdf_dir}")
        
        # Check if data has already been processed
        split_dir = os.path.join(SCRIPT_DIR, "../data/column_three/ligands_sdf_split")
        
        if os.path.exists(split_dir) and os.listdir(split_dir):
            # Count existing processed molecules
            existing_tranches = 0
            existing_molecules = 0
            for item in os.listdir(split_dir):
                item_path = os.path.join(split_dir, item)
                if os.path.isdir(item_path):
                    sdf_files = list(Path(item_path).glob("*.sdf"))
                    if sdf_files:
                        existing_tranches += 1
                        existing_molecules += len(sdf_files)
            
            logging.info(f"\n=== EXISTING PROCESSED DATA DETECTED ===")
            logging.info(f"✓ Found {existing_tranches} tranches with {existing_molecules:,} individual molecules")
            logging.info(f"📁 Location: {split_dir}")
            logging.info(f"\n🚀 SDF data is already processed and ready for docking!")
            
        elif successful_extractions > 0:
            # Split multi-molecule SDF files into individual molecules organized by tranche
            logging.info(f"\n--- Splitting extracted SDF files into individual molecules ---")
            total_molecules, failed_splits, tranche_count = split_sdf_files(
                sdf_dir, split_dir, max_workers=EXTRACTION_WORKERS)
            logging.info(f"\n=== SPLITTING SUMMARY ===")
            logging.info(f"✓ Total individual molecules created: {total_molecules}")
            logging.info(f"✗ Failed file splits: {failed_splits}")
            logging.info(f"📁 Tranches created: {tranche_count}")
            logging.info(f"📁 Individual SDF files organized by tranche: {split_dir}")
            
            if total_molecules > 0:
                logging.info(f"\n🚀 SUCCESS: {total_molecules} individual SDF molecules are ready for docking!")
            else:
                logging.warning("\n⚠️  WARNING: No molecules were successfully split.")
        else:
            logging.warning("\n⚠️  WARNING: No SDF files were successfully extracted. Check the download and extraction logs.")
                        
    except Exception as e:
        logging.error(f"Error during SDF data download: {e}")
        exit(1) 
