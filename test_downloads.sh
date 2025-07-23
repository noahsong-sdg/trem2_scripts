#!/bin/bash

echo "Testing different download methods from test.uri file"
echo "Working directory: $(pwd)"
echo "Date: $(date)"
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
echo "=== SUMMARY ==="
echo "Check the logs above to see which method worked:"
echo "- test_wget_3d.log"
echo "- test_curl_3d.log" 
echo "- test_mkdir_wget_3d.log"
echo ""
echo "Test completed at: $(date)" 
