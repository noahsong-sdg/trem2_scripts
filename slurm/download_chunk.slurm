#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --account=ucb-general
#SBATCH --job-name=download_split
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=download_split_%j.out
#SBATCH --error=download_split_%j.err

#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=noso3320@colorado.edu

SCRIPT_DIR="/scratch/alpine/noso3320/trem2/unidock/scripts"
DATA_DIR="/scratch/alpine/noso3320/trem2/unidock/data"
RESULTS_DIR="/scratch/alpine/noso3320/trem2/unidock/results"

module purge
module load gcc
module load cmake
module load anaconda

cd /scratch/alpine/noso3320/trem2/unidock
conda activate dock_env

echo "=== ZINC Download with Auto-Splitting ==="
echo "Started at: $(date)"
echo "Current directory: $(pwd)"

# Step 1: Check if rsync file exists
RSYNC_FILE="${DATA_DIR}/zinc22a.rsync"
if [ ! -f "${RSYNC_FILE}" ]; then
    echo "Error: Rsync file not found: ${RSYNC_FILE}"
    exit 1
fi

echo "Found rsync file: ${RSYNC_FILE}"
echo "File size: $(wc -l < ${RSYNC_FILE}) lines"

# Step 2: Split the rsync file into chunks
echo "Splitting rsync file into chunks..."
python scripts/split_rsync.py "${RSYNC_FILE}" 1000

if [ $? -ne 0 ]; then
    echo "Error: Failed to split rsync file"
    exit 1
fi

# Step 3: Count chunks and create array job script
CHUNK_DIR="${DATA_DIR}/zinc22a_chunks"
CHUNK_COUNT=$(ls ${CHUNK_DIR}/*.rsync | wc -l)
echo "Created ${CHUNK_COUNT} chunks"

# Step 4: Create array job script for processing chunks
ARRAY_SCRIPT="${SCRIPT_DIR}/slurm/process_chunks.slurm"
cat > "${ARRAY_SCRIPT}" << EOF
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --account=ucb-general
#SBATCH --job-name=process_chunk_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=process_chunk_%a_%j.out
#SBATCH --error=process_chunk_%a_%j.err
#SBATCH --array=1-${CHUNK_COUNT}

#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=noso3320@colorado.edu

SCRIPT_DIR="/scratch/alpine/noso3320/trem2/unidock/scripts"
DATA_DIR="/scratch/alpine/noso3320/trem2/unidock/data"

module purge
module load gcc
module load cmake
module load anaconda

cd /scratch/alpine/noso3320/trem2/unidock
conda activate dock_env

# Get the chunk file for this array job
CHUNK_DIR="\${DATA_DIR}/zinc22a_chunks"
CHUNK_FILES=(\$(ls \${CHUNK_DIR}/*.rsync 2>/dev/null | sort))
if [ \${#CHUNK_FILES[@]} -eq 0 ]; then
    echo "Error: No chunk files found in \${CHUNK_DIR}"
    exit 1
fi
CHUNK_FILE="\${CHUNK_FILES[\$((SLURM_ARRAY_TASK_ID-1))]}"
if [ -z "\${CHUNK_FILE}" ]; then
    echo "Error: No chunk file for task ID \${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Processing chunk: \${CHUNK_FILE}"
echo "Array task ID: \${SLURM_ARRAY_TASK_ID}"
echo "Started at: \$(date)"

# Run the download script with this chunk
python scripts/download.py "\${CHUNK_FILE}"

echo "Finished at: \$(date)"
echo "Exit code: \$?"
EOF

echo "Created array job script: ${ARRAY_SCRIPT}"

# Step 5: Submit the array job
echo "Submitting array job for ${CHUNK_COUNT} chunks..."
JOB_ID=$(sbatch "${ARRAY_SCRIPT}" | awk '{print $4}')
echo "Array job submitted with ID: ${JOB_ID}"

echo "=== Setup Complete ==="
echo "Split job finished at: $(date)"
echo "Array job ID: ${JOB_ID}"
echo "Monitor with: squeue -u \$USER"
echo "Check progress with: ls \${DATA_DIR}/zinc22a_chunks/"
