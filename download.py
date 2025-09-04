#!/usr/bin/env python3
"""
ZINC22A Ligand Download Script

This script parses rsync files and downloads ligand files from the ZINC22-3D database.
The script will look for rsync files in multiple locations for convenience.

Usage:
    python download.py                    # Uses zinc22a.rsync (default)
    python download.py myfile.rsync      # Uses specified rsync file
    python download.py --dry-run         # Test run without downloading

The script will automatically create the appropriate directory structure in ../data/
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import List, Dict, Optional
import time


class ZINC22ADownloader:
    """Handles downloading of ZINC22A ligand files using rsync commands."""
    
    def __init__(self, rsync_file: str, dry_run: bool = False):
        """
        Initialize the downloader.
        
        Args:
            rsync_file: Name or path of the rsync file to use
            dry_run: If True, only show what would be downloaded without actually downloading
        """
        self.script_dir = Path(__file__).parent
        self.dry_run = dry_run
        
        # Resolve rsync file path - check multiple locations
        self.rsync_file = self._resolve_rsync_file(rsync_file)
        
        # Auto-detect target directory based on rsync filename
        rsync_stem = Path(rsync_file).stem  # e.g., "zinc22a" from "zinc22a.rsync"
        self.target_dir = self.script_dir.parent / "data" / rsync_stem
        
        # Validate inputs
        if not self.rsync_file.exists():
            raise FileNotFoundError(f"Rsync file not found: {self.rsync_file}")
        
        # Create target directory if it doesn't exist
        self.target_dir.mkdir(parents=True, exist_ok=True)
    
    def _resolve_rsync_file(self, rsync_file: str) -> Path:
        """
        Resolve the rsync file path by checking multiple locations.
        
        Args:
            rsync_file: Name or path of the rsync file
            
        Returns:
            Resolved Path object to the rsync file
        """
        # If it's already a full path, use it directly
        if Path(rsync_file).is_absolute():
            return Path(rsync_file)
        
        # Check multiple locations in order of preference
        locations = [
            self.script_dir / rsync_file,                    # Current script directory
            self.script_dir.parent / "data" / rsync_file,   # ../data/ directory
            Path.cwd() / rsync_file,                        # Current working directory
        ]
        
        for location in locations:
            if location.exists():
                return location
        
        # If not found, return the first location for better error messages
        return locations[0]
        
    def parse_rsync_file(self) -> List[Dict[str, str]]:
        """
        Parse the rsync file and extract download commands.
        
        Returns:
            List of dictionaries containing parsed rsync command information
        """
        commands = []
        
        try:
            with open(self.rsync_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    
                    # Skip empty lines and comments
                    if not line or line.startswith('#'):
                        continue
                    
                    # Parse rsync commands
                    if line.startswith('rsync '):
                        command_info = self._parse_rsync_line(line, line_num)
                        if command_info:
                            commands.append(command_info)
                    
                    # Parse mkdir and pushd commands for context
                    elif line.startswith('mkdir ') or line.startswith('pushd '):
                        commands.append({
                            'type': 'command',
                            'command': line,
                            'line_num': line_num,
                            'description': f'Directory operation: {line}'
                        })
        
        except Exception as e:
            print(f"Error parsing rsync file: {e}")
            raise
        
        print(f"Parsed {len(commands)} commands from {self.rsync_file.name}")
        return commands
    
    def _parse_rsync_line(self, line: str, line_num: int) -> Optional[Dict[str, str]]:
        """
        Parse a single rsync command line.
        
        Args:
            line: The rsync command line
            line_num: Line number for error reporting
            
        Returns:
            Dictionary with parsed command information or None if parsing fails
        """
        try:
            # Extract key components from rsync command
            parts = line.split()
            
            # Find the source (rsync://...)
            source = None
            for part in parts:
                if part.startswith('rsync://'):
                    source = part
                    break
            
            if not source:
                print(f"Warning: Could not parse source from line {line_num}")
                return None
            
            # Extract include/exclude patterns
            includes = []
            excludes = []
            for i, part in enumerate(parts):
                if part == '--include=':
                    if i + 1 < len(parts):
                        includes.append(parts[i + 1])
                elif part.startswith('--include='):
                    includes.append(part.split('=', 1)[1])
                elif part == '--exclude=':
                    if i + 1 < len(parts):
                        excludes.append(parts[i + 1])
                elif part.startswith('--exclude='):
                    excludes.append(part.split('=', 1)[1])
            
            return {
                'type': 'rsync',
                'command': line,
                'line_num': line_num,
                'source': source,
                'includes': includes,
                'excludes': excludes,
                'description': f'Download from {source} with {len(includes)} includes, {len(excludes)} excludes'
            }
            
        except Exception as e:
            print(f"Error parsing line {line_num}: {e}")
            return None
    
    def execute_commands(self, commands: List[Dict[str, str]]) -> None:
        """
        Execute the parsed commands.
        
        Args:
            commands: List of parsed command dictionaries
        """
        total_commands = len([c for c in commands if c['type'] == 'rsync'])
        completed = 0
        failed = 0
        
        print(f"Starting execution of {total_commands} download commands...")
        
        # Change to target directory
        original_cwd = os.getcwd()
        os.chdir(self.target_dir)
        
        try:
            for command_info in commands:
                if command_info['type'] == 'command':
                    # Handle mkdir, pushd, etc.
                    self._execute_system_command(command_info)
                elif command_info['type'] == 'rsync':
                    # Handle rsync downloads
                    success = self._execute_rsync_command(command_info)
                    if success:
                        completed += 1
                        print(f"✓ Downloaded ({completed}/{total_commands}) from {command_info['source']}")
                    else:
                        failed += 1
                        print(f"✗ Failed ({failed}/{total_commands}) from {command_info['source']}")
                    
                    # Progress update every 10 commands
                    if (completed + failed) % 10 == 0:
                        print(f"Progress: {completed + failed}/{total_commands} completed")
        
        finally:
            # Restore original working directory
            os.chdir(original_cwd)
        
        print(f"\n=== DOWNLOAD SUMMARY ===")
        print(f"✓ Successful downloads: {completed}")
        print(f"✗ Failed downloads: {failed}")
        print(f"📁 Files saved to: {self.target_dir}")
    
    def _execute_system_command(self, command_info: Dict[str, str]) -> None:
        """Execute a system command (mkdir, pushd, etc.)."""
        command = command_info['command']
        line_num = command_info['line_num']
        
        if self.dry_run:
            print(f"[DRY RUN] Would execute: {command}")
            return
        
        try:
            if command.startswith('mkdir '):
                # Extract directory name and create it
                dir_name = command.split(' ', 1)[1]
                Path(dir_name).mkdir(parents=True, exist_ok=True)
                print(f"Created directory: {dir_name}")
            elif command.startswith('pushd '):
                # Extract directory name and change to it
                dir_name = command.split(' ', 1)[1]
                os.chdir(dir_name)
                print(f"Changed to directory: {dir_name}")
            else:
                print(f"Unknown command type at line {line_num}: {command}")
        
        except Exception as e:
            print(f"Error executing command at line {line_num}: {e}")
    
    def _execute_rsync_command(self, command_info: Dict[str, str]) -> bool:
        """
        Execute a rsync download command.
        
        Args:
            command_info: Parsed rsync command information
            
        Returns:
            True if successful, False otherwise
        """
        command = command_info['command']
        line_num = command_info['line_num']
        source = command_info['source']
        
        if self.dry_run:
            print(f"[DRY RUN] Would download from: {source}")
            return True
        
        try:
            print(f"Downloading from {source}...")
            
            # Execute rsync command
            result = subprocess.run(
                command.split(),
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )
            
            if result.returncode == 0:
                return True
            else:
                print(f"Rsync failed for {source}: {result.stderr}")
                return False
        
        except subprocess.TimeoutExpired:
            print(f"Rsync timed out for {source}")
            return False
        except Exception as e:
            print(f"Error executing rsync for {source}: {e}")
            return False
    
    def run(self) -> None:
        """Main execution method."""
        start_time = time.time()
        
        try:
            print(f"=== ZINC22A LIGAND DOWNLOAD SCRIPT ===")
            print(f"Rsync file: {self.rsync_file}")
            print(f"Target directory: {self.target_dir}")
            print(f"Dry run: {self.dry_run}")
            print()
            
            # Parse the rsync file
            commands = self.parse_rsync_file()
            
            if not commands:
                print("No commands found in rsync file")
                return
            
            # Execute the commands
            self.execute_commands(commands)
            
            elapsed_time = time.time() - start_time
            print(f"\nDownload process completed in {elapsed_time:.2f} seconds")
            
        except Exception as e:
            print(f"Download process failed: {e}")
            raise


def main():
    """Main entry point."""
    # Simple argument parsing - just check for dry-run flag
    dry_run = '--dry-run' in sys.argv
    
    # Get rsync filename from command line or use default
    rsync_file = 'zinc22a.rsync'  # default
    for arg in sys.argv[1:]:
        if not arg.startswith('--') and arg.endswith('.rsync'):
            rsync_file = arg
            break
    
    try:
        # Create and run the downloader
        downloader = ZINC22ADownloader(
            rsync_file=rsync_file,
            dry_run=dry_run
        )
        
        downloader.run()
        
    except KeyboardInterrupt:
        print("\nDownload interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
