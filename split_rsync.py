#!/usr/bin/env python3
"""
Rsync File Splitter

This script splits a large rsync file into smaller chunks for easier processing.
Each chunk will be a complete, valid rsync file that can be processed independently.

Usage:
    python split_rsync.py zinc22a.rsync 1000    # Split into chunks of 1000 lines each
    python split_rsync.py zinc22a.rsync 5000    # Split into chunks of 5000 lines each
"""

import os
import sys
from pathlib import Path


def split_rsync_file(rsync_file, lines_per_chunk=1000):
    """
    Split a large rsync file into smaller chunks.
    
    Args:
        rsync_file (str): Path to the rsync file to split
        lines_per_chunk (int): Number of rsync commands per chunk
    """
    rsync_path = Path(rsync_file)
    
    if not rsync_path.exists():
        print(f"Error: Rsync file not found: {rsync_file}")
        return False
    
    # Create output directory for chunks
    output_dir = rsync_path.parent / f"{rsync_path.stem}_chunks"
    output_dir.mkdir(exist_ok=True)
    
    print(f"Splitting {rsync_file} into chunks of {lines_per_chunk} rsync commands each...")
    print(f"Output directory: {output_dir}")
    
    chunk_num = 1
    current_chunk_lines = 0
    current_chunk_file = None
    
    try:
        with open(rsync_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                
                # Start new chunk if needed
                if current_chunk_lines == 0:
                    if current_chunk_file:
                        current_chunk_file.close()
                    
                    chunk_filename = output_dir / f"{rsync_path.stem}_chunk_{chunk_num:03d}.rsync"
                    current_chunk_file = open(chunk_filename, 'w')
                    print(f"Creating chunk {chunk_num}: {chunk_filename.name}")
                
                # Write the line to current chunk
                current_chunk_file.write(line + '\n')
                current_chunk_lines += 1
                
                # Check if we need to start a new chunk
                if line.startswith('rsync ') and current_chunk_lines >= lines_per_chunk:
                    current_chunk_file.close()
                    print(f"  Chunk {chunk_num} complete: {current_chunk_lines} rsync commands")
                    chunk_num += 1
                    current_chunk_lines = 0
                    current_chunk_file = None
        
        # Close the last chunk if it exists
        if current_chunk_file:
            current_chunk_file.close()
            print(f"  Chunk {chunk_num} complete: {current_chunk_lines} rsync commands")
        
        print(f"\nSplit complete! Created {chunk_num} chunks in {output_dir}")
        print(f"Each chunk can be processed independently with download.py")
        
        return True
        
    except Exception as e:
        print(f"Error splitting rsync file: {e}")
        if current_chunk_file:
            current_chunk_file.close()
        return False


def main():
    if len(sys.argv) < 2:
        print("Usage: python split_rsync.py <rsync_file> [lines_per_chunk]")
        print("Example: python split_rsync.py zinc22a.rsync 1000")
        sys.exit(1)
    
    rsync_file = sys.argv[1]
    lines_per_chunk = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    
    print(f"Rsync file: {rsync_file}")
    print(f"Lines per chunk: {lines_per_chunk}")
    
    success = split_rsync_file(rsync_file, lines_per_chunk)
    
    if success:
        print("\nNext steps:")
        print("1. Check the chunks in the _chunks directory")
        print("2. Process each chunk with: python download.py chunk_name.rsync")
        print("3. Or create separate SLURM jobs for each chunk")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
