#!/bin/bash
#SBATCH --job-name=test_zinc_downloads
#SBATCH --partition=amilan
#SBATCH --account=ucb-general
#SBATCH --qos=normal
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=test_zinc_downloads_%j.out
#SBATCH --error=test_zinc_downloads_%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=noso3320@colorado.edu

# Load necessary modules
module purge
module load gcc

# Set working directory
WORK_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock"
cd $WORK_DIR

echo "Testing different download methods from test.uri file"
echo "Working directory: $(pwd)"
echo "Date: $(date)"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo ""

# Clean up any previous test files
rm -rf AC/
rm -f ACAAML.xaa.pdbqt.gz
rm -f test_*.log

# Test 1: Direct wget with 3D path
echo "=== TEST 1: Direct wget with /3D/ path ==="
echo "Command: wget http://files.docking.org/3D/AC/AAML/ACAAML.xaa.pdbqt.gz"
if wget --timeout=120 --connect-timeout=30 --tries=2 http://files.docking.org/3D/AC/AAML/ACAAML.xaa.pdbqt.gz 2>&1 | tee test_wget_3d.log; then
    echo "✓ SUCCESS: Direct wget with /3D/ worked"
    ls -la ACAAML.xaa.pdbqt.gz
    rm -f ACAAML.xaa.pdbqt.gz
else
    echo "✗ FAILED: Direct wget with /3D/"
fi
echo ""

# Test 2: curl with create-dirs
echo "=== TEST 2: curl with --create-dirs ==="
echo "Command: curl --remote-time --fail --create-dirs -o AC/AAML/ACAAML.xaa.pdbqt.gz http://files.docking.org/3D/AC/AAML/ACAAML.xaa.pdbqt.gz"
if curl --connect-timeout 30 --max-time 120 --remote-time --fail --create-dirs -o AC/AAML/ACAAML.xaa.pdbqt.gz http://files.docking.org/3D/AC/AAML/ACAAML.xaa.pdbqt.gz 2>&1 | tee test_curl_3d.log; then
    echo "✓ SUCCESS: curl with --create-dirs worked"
    ls -la AC/AAML/ACAAML.xaa.pdbqt.gz
    rm -rf AC/
else
    echo "✗ FAILED: curl with --create-dirs"
fi
echo ""

# Test 3: mkdir + wget
echo "=== TEST 3: mkdir + wget ==="
echo "Command: mkdir -pv AC/AAML && wget http://files.docking.org/3D/AC/AAML/ACAAML.xaa.pdbqt.gz -O AC/AAML/ACAAML.xaa.pdbqt.gz"
if mkdir -pv AC/AAML && wget --timeout=120 --connect-timeout=30 --tries=2 http://files.docking.org/3D/AC/AAML/ACAAML.xaa.pdbqt.gz -O AC/AAML/ACAAML.xaa.pdbqt.gz 2>&1 | tee test_mkdir_wget_3d.log; then
    echo "✓ SUCCESS: mkdir + wget worked"
    ls -la AC/AAML/ACAAML.xaa.pdbqt.gz
    rm -rf AC/
else
    echo "✗ FAILED: mkdir + wget"
fi
echo ""

# Test 4: Test connectivity to both 2D and 3D paths
echo "=== TEST 4: Compare /2D/ vs /3D/ paths ==="
echo "Testing HEAD requests..."

echo "Testing /2D/ path:"
curl -I --connect-timeout 30 --max-time 60 http://files.docking.org/2D/AC/AAML/ACAAML.xaa.pdbqt.gz 2>&1 | head -3

echo "Testing /3D/ path:"
curl -I --connect-timeout 30 --max-time 60 http://files.docking.org/3D/AC/AAML/ACAAML.xaa.pdbqt.gz 2>&1 | head -3

echo ""

# Test 5: Test our database index path structure
echo "=== TEST 5: Database index path structure ==="
echo "Our database index has paths like: AB/AAML/ABAAML.xaa.pdbqt.gz"
echo "Testing if this works with /3D/ base:"

echo "Testing /3D/ + database index path:"
curl -I --connect-timeout 30 --max-time 60 http://files.docking.org/3D/AB/AAML/ABAAML.xaa.pdbqt.gz 2>&1 | head -3

echo "Testing /2D/ + database index path:"
curl -I --connect-timeout 30 --max-time 60 http://files.docking.org/2D/AB/AAML/ABAAML.xaa.pdbqt.gz 2>&1 | head -3

echo ""

# Test 6: Test cluster network connectivity
echo "=== TEST 6: Network connectivity tests ==="
echo "Testing ping to files.docking.org:"
ping -c 3 files.docking.org

echo "Testing basic HTTP connectivity:"
curl -I --connect-timeout 30 http://files.docking.org/ 2>&1 | head -3

echo ""
echo "=== SUMMARY ==="
echo "Check the logs above to see which method worked:"
echo "- test_wget_3d.log"
echo "- test_curl_3d.log" 
echo "- test_mkdir_wget_3d.log"
echo ""
echo "Key findings:"
echo "1. Which download method worked: [see above]"
echo "2. /2D/ vs /3D/ path: [see TEST 4 results]"
echo "3. Database index path compatibility: [see TEST 5 results]"
echo "4. Network connectivity: [see TEST 6 results]"
echo ""
echo "Test completed at: $(date)"

exit 0 
