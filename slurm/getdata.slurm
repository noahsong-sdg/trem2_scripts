#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --account=ucb-general
#SBATCH --job-name=get_column_one
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --output=c1_%j.out
#SBATCH --error=c1_%j.err

#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=noso3320@colorado.edu


SCRIPT_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock/trem2_scripts"
DATA_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock/data"
RESULTS_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock/results"

module purge
module load gcc
module load cmake
module load anaconda

cd /scratch/alpine/noso3320/trem2/trem2_unidock
conda activate dock_env

python trem2_scripts/getdata_sdf.py
