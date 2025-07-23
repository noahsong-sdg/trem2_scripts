#!/usr/bin/env python3
import requests
import os
import gzip
import time
import random
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configuration
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_LIGANDS_DIR = os.path.join(ROOT_DIR, "../data/test/ligands_raw")
URI_FILE = os.path.join(ROOT_DIR, "test.uri")

def download_zinc_subset(url, output_dir, filename=None):
    """Download a single ZINC file."""
    os.makedirs(output_dir, exist_ok=True)
    
    if filename is None:
        filename = url.split('/')[-1]
    
    filepath = os.path.join(output_dir, filename)
    time.sleep(random.uniform(0.5, 2.0))  # Rate limiting
    
    print(f"Downloading {filename} from {url}")
    start_time = time.time()
    
    # Longer timeouts for large files and potentially slow connections
    response = requests.get(url, stream=True, timeout=(30, 600),  # (connect, read) timeout
                          headers={'User-Agent': 'HTS-Pipeline/1.0'})
    response.raise_for_status()
    
    # Get file size if available
    total_size = int(response.headers.get('content-length', 0))
    downloaded = 0
    
    with open(filepath, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
                downloaded += len(chunk)
    
    elapsed = time.time() - start_time
    speed = downloaded / elapsed / 1024 if elapsed > 0 else 0  # KB/s
    print(f"Downloaded: {filename} ({downloaded:,} bytes in {elapsed:.1f}s, {speed:.1f} KB/s)")
    return filepath

def download_all_from_uri_file(uri_file_path, base_output_dir, max_workers=3):
    """Download all URLs from URI file using parallel processing."""
    if not os.path.exists(uri_file_path):
        print(f"Error: URI file not found at {uri_file_path}")
        return 0, 0

    # Read URLs
    with open(uri_file_path, 'r') as f:
        urls = [line.strip() for line in f if line.strip() and not line.startswith("#")]
    
    if not urls:
        print("No valid URLs found")
        return 0, 0
    
    print(f"Downloading {len(urls)} files with {max_workers} workers...")
    
    successful = 0
    failed = 0
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit downloads directly to download_zinc_subset
        futures = {executor.submit(download_zinc_subset, url, base_output_dir): url for url in urls}
        
        for future in as_completed(futures):
            try:
                if future.result():
                    successful += 1
                else:
                    failed += 1
            except Exception as e:
                failed += 1
                print(f"Download failed: {e}")
    
    print(f"Complete: {successful} successful, {failed} failed")
    return successful, failed

def extract_pdbqt_files(raw_dir, output_dir):
    """Extract .pdbqt.gz files to .pdbqt files."""
    os.makedirs(output_dir, exist_ok=True)
    
    gz_files = list(Path(raw_dir).glob("*.pdbqt.gz"))
    if not gz_files:
        print(f"No .pdbqt.gz files found in {raw_dir}")
        return 0, 0
    
    print(f"Extracting {len(gz_files)} files...")
    
    successful = 0
    for gz_file in gz_files:
        try:
            output_file = Path(output_dir) / gz_file.stem
            with gzip.open(gz_file, 'rb') as f_in:
                with open(output_file, 'wb') as f_out:
                    f_out.write(f_in.read())
            successful += 1
        except Exception as e:
            print(f"Extraction failed for {gz_file.name}: {e}")
    
    print(f"Extracted {successful}/{len(gz_files)} files")
    return successful, len(gz_files) - successful

def split_pdbqt_files(input_dir, output_dir):
    """Split multi-molecule PDBQT files into individual files by tranche."""
    os.makedirs(output_dir, exist_ok=True)
    
    pdbqt_files = list(Path(input_dir).glob("*.pdbqt"))
    if not pdbqt_files:
        print(f"No .pdbqt files found in {input_dir}")
        return 0
    
    print(f"Splitting {len(pdbqt_files)} files...")
    
    total_molecules = 0
    for pdbqt_file in pdbqt_files:
        tranche_name = pdbqt_file.stem
        tranche_dir = Path(output_dir) / tranche_name
        tranche_dir.mkdir(exist_ok=True)
        
        molecule_count = 0
        current_molecule = []
        
        with open(pdbqt_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                if line.startswith('MODEL'):
                    if current_molecule:
                        _save_molecule(current_molecule, tranche_dir, molecule_count)
                        molecule_count += 1
                    current_molecule = []
                elif line.startswith('ENDMDL'):
                    if current_molecule:
                        _save_molecule(current_molecule, tranche_dir, molecule_count)
                        molecule_count += 1
                        current_molecule = []
                else:
                    if not (line.startswith('MODEL') or line.startswith('ENDMDL')):
                        current_molecule.append(line)
        
        if current_molecule:
            _save_molecule(current_molecule, tranche_dir, molecule_count)
            molecule_count += 1
            
        total_molecules += molecule_count
        print(f"Split {pdbqt_file.name}: {molecule_count} molecules")
    
    print(f"Total molecules: {total_molecules}")
    return total_molecules

def _save_molecule(molecule_lines, tranche_dir, molecule_index):
    """Save individual molecule to PDBQT file."""
    if not molecule_lines:
        return
        
    filename = f"molecule_{molecule_index:06d}.pdbqt"
    output_path = tranche_dir / filename
    
    with open(output_path, 'w') as f:
        for line in molecule_lines:
            f.write(line + '\n')

if __name__ == "__main__":
    print("=== ZINC Download Pipeline ===")

    successful_downloads, failed_downloads = download_all_from_uri_file(
        URI_FILE, RAW_LIGANDS_DIR, max_workers=4)
    
    if successful_downloads == 0:
        print("No files downloaded successfully")
        exit(1)
    
    pdbqt_dir = os.path.join(ROOT_DIR, "../data/test/ligands_pdbqt")
    extract_pdbqt_files(RAW_LIGANDS_DIR, pdbqt_dir)
    
    # Split
    split_dir = os.path.join(ROOT_DIR, "../data/test/ligands_pdbqt_split")
    total_molecules = split_pdbqt_files(pdbqt_dir, split_dir)
    
    print(f"Pipeline complete: {total_molecules} molecules ready at {split_dir}")
