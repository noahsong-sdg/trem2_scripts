#!/usr/bin/env python3
"""
Debug script to check download status and diagnose issues.
"""

import os
import sys
from pathlib import Path


def check_download_status():
    """Check the current status of downloads."""
    print("=== DOWNLOAD STATUS CHECK ===")
    
    # Check data directory structure
    data_dir = Path("../data")
    print(f"Data directory: {data_dir.absolute()}")
    
    if not data_dir.exists():
        print("❌ Data directory does not exist!")
        return
    
    # Check for rsync file
    rsync_file = data_dir / "zinc22a.rsync"
    if rsync_file.exists():
        print(f"✅ Found rsync file: {rsync_file}")
        print(f"   Size: {rsync_file.stat().st_size / (1024*1024):.1f} MB")
    else:
        print(f"❌ Rsync file not found: {rsync_file}")
    
    # Check for chunks directory
    chunks_dir = data_dir / "zinc22a_chunks"
    if chunks_dir.exists():
        chunk_files = list(chunks_dir.glob("*.rsync"))
        print(f"✅ Found chunks directory: {chunks_dir}")
        print(f"   Number of chunk files: {len(chunk_files)}")
        
        if chunk_files:
            # Check first few chunks
            for i, chunk_file in enumerate(sorted(chunk_files)[:3]):
                print(f"   Chunk {i+1}: {chunk_file.name} ({chunk_file.stat().st_size} bytes)")
                
                # Count rsync commands in first chunk
                with open(chunk_file, 'r') as f:
                    rsync_count = sum(1 for line in f if line.strip().startswith('rsync '))
                print(f"     Rsync commands: {rsync_count}")
    else:
        print(f"❌ Chunks directory not found: {chunks_dir}")
    
    # Check for download outputs
    print("\n=== DOWNLOAD OUTPUTS ===")
    
    # Look for any downloaded data
    for subdir in ["zinc22a"]:
        target_dir = data_dir / subdir
        if target_dir.exists():
            print(f"✅ Found target directory: {target_dir}")
            
            # Check subdirectories
            for subsubdir in ["ligands_raw", "ligands_sdf", "ligands_sdf_split"]:
                check_dir = target_dir / subsubdir
                if check_dir.exists():
                    files = list(check_dir.iterdir())
                    print(f"   {subsubdir}: {len(files)} files")
                    if files:
                        print(f"     Example: {files[0].name}")
                else:
                    print(f"   {subsubdir}: Not found")
        else:
            print(f"❌ Target directory not found: {target_dir}")
    
    # Check SLURM output files
    print("\n=== SLURM OUTPUTS ===")
    current_dir = Path(".")
    slurm_files = list(current_dir.glob("*download*.out")) + list(current_dir.glob("*download*.err"))
    
    if slurm_files:
        print(f"Found {len(slurm_files)} SLURM output files:")
        for slurm_file in sorted(slurm_files):
            print(f"   {slurm_file.name}")
            if slurm_file.stat().st_size > 0:
                print(f"     Size: {slurm_file.stat().st_size} bytes")
                # Show last few lines
                try:
                    with open(slurm_file, 'r') as f:
                        lines = f.readlines()
                        if lines:
                            print(f"     Last line: {lines[-1].strip()}")
                except:
                    print(f"     Could not read file")
            else:
                print(f"     Empty file")
    else:
        print("❌ No SLURM output files found")


def test_chunk_processing():
    """Test processing a single chunk."""
    print("\n=== TESTING CHUNK PROCESSING ===")
    
    chunks_dir = Path("../data/zinc22a_chunks")
    if not chunks_dir.exists():
        print("❌ No chunks directory found")
        return
    
    chunk_files = list(chunks_dir.glob("*.rsync"))
    if not chunk_files:
        print("❌ No chunk files found")
        return
    
    # Test with first chunk
    test_chunk = sorted(chunk_files)[0]
    print(f"Testing with chunk: {test_chunk.name}")
    
    # Check chunk content
    with open(test_chunk, 'r') as f:
        lines = f.readlines()
        rsync_lines = [line for line in lines if line.strip().startswith('rsync ')]
        print(f"   Total lines: {len(lines)}")
        print(f"   Rsync commands: {len(rsync_lines)}")
        
        if rsync_lines:
            print(f"   First rsync command: {rsync_lines[0].strip()[:100]}...")
    
    # Test download script with dry run
    print(f"\nTesting download script with dry run...")
    import subprocess
    
    try:
        result = subprocess.run([
            "python", "download.py", str(test_chunk), "--dry-run"
        ], capture_output=True, text=True, timeout=60)
        
        print(f"   Return code: {result.returncode}")
        if result.stdout:
            print(f"   STDOUT: {result.stdout[:500]}...")
        if result.stderr:
            print(f"   STDERR: {result.stderr[:500]}...")
            
    except subprocess.TimeoutExpired:
        print("   ❌ Command timed out (60 seconds)")
    except Exception as e:
        print(f"   ❌ Error: {e}")


if __name__ == "__main__":
    check_download_status()
    test_chunk_processing()
