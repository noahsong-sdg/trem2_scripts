#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --qos=preemptable
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24 #16 ran out of mem
#SBATCH --job-name=mcdock_b3
#SBATCH --output=mcdockb3_%j.out
#SBATCH --error=mcdockb3_%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=noso3320@colorado.edu

module purge
module load gcc
module load cuda/11.8
module load cmake
module load anaconda

SCRIPT_DIR="/scratch/alpine/noso3320/trem2/unidock/scripts"
cd /scratch/alpine/noso3320/trem2/unidock
conda activate dock_env
python $SCRIPT_DIR/mcdockv2_b3.py
