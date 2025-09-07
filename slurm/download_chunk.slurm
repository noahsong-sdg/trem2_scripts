#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --account=ucb-general
#SBATCH --job-name=download_chunk_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=download_chunk_%a_%j.out
#SBATCH --error=download_chunk_%a_%j.err
#SBATCH --array=1-10  # Adjust based on number of chunks

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

# Get the chunk file for this array job
CHUNK_DIR="${DATA_DIR}/zinc22a_chunks"
CHUNK_FILES=($(ls ${CHUNK_DIR}/*.rsync | sort))
CHUNK_FILE="${CHUNK_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

echo "Processing chunk: ${CHUNK_FILE}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Started at: $(date)"

# Run the download script with this chunk
python scripts/download.py "${CHUNK_FILE}"

echo "Finished at: $(date)"
echo "Exit code: $?"
