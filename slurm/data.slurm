#!/bin/bash
#SBATCH --job-name=zinc_download
#SBATCH --partition=amilan
#SBATCH --account=ucb-general
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=zinc_download_%j.out
#SBATCH --error=zinc_download_%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=noso3320@colorado.edu

# Load necessary modules
module purge
module load gcc
module load cmake

# Set working directory paths to match project structure
SCRIPT_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock/scripts"
DATA_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock/data"
WORK_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock"
cd $WORK_DIR

# Database index file (relative paths to ZINC files)
DB_INDEX_FILE="$WORK_DIR/c1_dbindex.database_index"

# Try different possible base URLs for ZINC
ZINC_BASE_URLS=(
    "http://files.docking.org/2D"
    "http://files.docking.org"
    "https://files.docking.org/2D"
    "https://files.docking.org"
    "http://zinc15.docking.org/2D"
    "https://zinc15.docking.org/2D"
)

# Check if database index file exists
if [ ! -f "$DB_INDEX_FILE" ]; then
    echo "ERROR: $DB_INDEX_FILE not found! Please ensure the file is in the scripts directory."
    exit 1
fi

# Test connectivity to files.docking.org
echo "Testing connectivity to files.docking.org..."
if ping -c 3 files.docking.org > /dev/null 2>&1; then
    echo "✓ Can ping files.docking.org"
else
    echo "⚠️  Cannot ping files.docking.org - may still work via HTTP"
fi

# Test basic HTTP connectivity
echo "Testing HTTP connectivity..."
if wget --spider --timeout=30 http://files.docking.org > /dev/null 2>&1; then
    echo "✓ Can reach files.docking.org via HTTP"
else
    echo "⚠️  Cannot reach files.docking.org via HTTP - downloads may fail"
fi

# Function to download a file with retries
download_with_retry() {
    local file_path="$1"
    local filename=$(basename "$file_path")
    
    echo "Trying to download: $filename"
    echo "File path in index: $file_path"
    
    # Try each base URL until one works
    for base_url in "${ZINC_BASE_URLS[@]}"; do
        local url="${base_url}/${file_path}"
        echo "  Testing base URL: $base_url"
        echo "  Full URL: $url"
        
        # First, test if the URL exists with HEAD request
        echo "  Testing URL accessibility..."
        local head_response=$(curl -I --connect-timeout 30 --max-time 60 "$url" 2>&1 | head -1)
        echo "  HEAD response: $head_response"
        
        if echo "$head_response" | grep -q "200 OK"; then
            echo "  ✓ URL accessible, proceeding with download..."
            
            # Try downloading with this base URL
            local max_attempts=2
            local attempt=1
            
            while [ $attempt -le $max_attempts ]; do
                echo "    Attempt $attempt/$max_attempts for $filename"
                
                # Try wget first
                if wget --continue --timeout=300 --connect-timeout=30 --read-timeout=180 --tries=1 --progress=bar:force "$url" 2>/dev/null; then
                    echo "    ✓ Successfully downloaded $filename with wget from $base_url"
                    return 0
                fi
                
                # Try curl as fallback
                if curl -L -C - --connect-timeout 30 --max-time 300 --retry 1 -o "$filename" "$url" 2>/dev/null; then
                    echo "    ✓ Successfully downloaded $filename with curl from $base_url"
                    return 0
                fi
                
                attempt=$((attempt + 1))
                if [ $attempt -le $max_attempts ]; then
                    echo "    Retrying in 30s..."
                    sleep 30
                fi
            done
            
            echo "  ✗ Download failed with working URL: $url"
        else
            echo "  ✗ URL not accessible: $url"
        fi
         done
    
    echo "✗ Failed to download $filename - no working base URL found"
    return 1
}

# Count total files in database index
TOTAL_FILES=$(wc -l < "$DB_INDEX_FILE")
echo "Found $TOTAL_FILES files to download from $DB_INDEX_FILE"

# TESTING: Limit to first 3 files for quick test
TEST_LIMIT=3
echo "🧪 TESTING MODE: Processing only first $TEST_LIMIT files (out of $TOTAL_FILES total)"

# Download each file from database index
FAILED_DOWNLOADS=()
CURRENT_COUNT=0

while IFS= read -r file_path; do
    # Skip empty lines
    if [ -z "$file_path" ]; then
        continue
    fi
    
    CURRENT_COUNT=$((CURRENT_COUNT + 1))
    echo "[$CURRENT_COUNT/$TEST_LIMIT] Processing: $file_path"
    
    if ! download_with_retry "$file_path"; then
        FAILED_DOWNLOADS+=("$file_path")
    fi
    
    # Break after test limit
    if [ $CURRENT_COUNT -ge $TEST_LIMIT ]; then
        echo "🧪 TESTING: Stopping after $TEST_LIMIT files"
        break
    fi
    
done < "$DB_INDEX_FILE"

# Copy downloaded files to data directory
echo "Copying downloaded files to data directory..."
mkdir -p "$DATA_DIR"
for file in *.pdbqt.gz; do
    if [ -f "$file" ]; then
        cp "$file" "$DATA_DIR/"
        echo "Copied: $file"
    fi
done

# Final report
echo ""
echo "========== DOWNLOAD SUMMARY =========="
echo "Total files requested: $CURRENT_COUNT"
echo "Successfully downloaded: $((CURRENT_COUNT - ${#FAILED_DOWNLOADS[@]}))"
echo "Failed downloads: ${#FAILED_DOWNLOADS[@]}"

if [ ${#FAILED_DOWNLOADS[@]} -gt 0 ]; then
    echo ""
    echo "Failed files:"
    for failed in "${FAILED_DOWNLOADS[@]}"; do
        echo "  - $failed"
    done
fi

echo ""
echo "Download job completed at: $(date)"
echo "Files are available in: $DATA_DIR"

exit 0
