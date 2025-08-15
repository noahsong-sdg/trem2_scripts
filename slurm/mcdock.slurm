#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=aa100 # aa100
#SBATCH --account=ucb-general
#SBATCH --job-name=mcdock

#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=16
#SBATCH --output=mcdock_%j.out
#SBATCH --error=mcdock_%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=noso3320@colorado.edu


SCRIPT_DIR="/scratch/alpine/noso3320/trem2/unidock/scripts"
DATA_DIR="/scratch/alpine/noso3320/trem2/unidock/data/column_one"
RESULTS_DIR="/scratch/alpine/noso3320/trem2/unidock/results"

module purge
module load gcc
module load cuda/11.8
module load cmake
module load anaconda

cd /scratch/alpine/noso3320/trem2/unidock
conda activate dock_env
echo "activated env, running script"
python $SCRIPT_DIR/mcdockv2.py
