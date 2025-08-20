#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --account=ucb-general
#SBATCH --job-name=get_column_fouri
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=c4i_%j.out
#SBATCH --error=c4i_%j.err

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

python scripts/getdata_sdf.py
