#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=aa100 
#SBATCH --account=ucb-general
#SBATCH --job-name=unidockdock

#SBATCH --nodes=1
#SBATCH --ntasks=12

#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=4      

#SBATCH --output=unidock_hts_%j.out
#SBATCH --error=unidock_hts_%j.err

#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=noso3320@colorado.edu

SCRIPT_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock/scripts"
DATA_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock/data"
RESULTS_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock/results"

module purge
module load gcc
module load cuda/11.8
module load cmake
module load anaconda

echo "=== UniDock High-Throughput Screening Workflow ==="
cd /scratch/alpine/noso3320/trem2/trem2_unidock
conda activate /projects/noso3320/hts_env

# --- IMPORTANT: Set the array size ---
# Run `wc -l ligand_manifest.txt` to find the total number of ligands.
# If it's 450,000, you would set the array to 1-450000.
#SBATCH --array=1-450000

# --- Define master paths ---
# Using environment variables makes the script clean and easy to modify.
export SCRIPT_DIR="$(pwd)"
export BASE_DIR="${SCRIPT_DIR}/.."
export RECEPTOR_FILE="${BASE_DIR}/data/receptor/cluster1_receptor.pdbqt"
export LIGAND_MANIFEST="${SCRIPT_DIR}/ligand_manifest.txt"
export MAIN_OUTPUT_DIR="${BASE_DIR}/results/array_outputs"
export LOG_DIR="${BASE_DIR}/results/array_logs"

# Create output directories if they don't exist
mkdir -p "${MAIN_OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"


# --- Core Job Array Logic ---
# Get the full path to the ligand for this specific job task.
# SLURM_ARRAY_TASK_ID is a special variable set by Slurm (e.g., 1, 2, 3...)
# 'sed' reads our manifest and picks the Nth line, where N is the task ID.
LIGAND_FILE_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${LIGAND_MANIFEST}")

# Check if the ligand file was found.
if [ -z "${LIGAND_FILE_PATH}" ]; then
    echo "Error: Could not find ligand for task ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Call our new, simple Python script with all parameters
python3 "${SCRIPT_DIR}/run_single_dock.py" \
    --receptor_path "${RECEPTOR_FILE}" \
    --ligand_path "${LIGAND_FILE_PATH}" \
    --output_dir "${MAIN_OUTPUT_DIR}" \
    --center_x 42.328 \
    --center_y 28.604 \
    --center_z 21.648 \
    --size_x 30.0 \
    --size_y 30.0 \
    --size_z 30.0

EXIT_CODE=$?

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "Job Task ${SLURM_ARRAY_TASK_ID} completed successfully."
else
    echo "Job Task ${SLURM_ARRAY_TASK_ID} failed with exit code ${EXIT_CODE}."
fi

exit ${EXIT_CODE}

    