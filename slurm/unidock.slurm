#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=aa100 # aa100
#SBATCH --account=ucb-general
#SBATCH --job-name=unidockdock

#SBATCH --nodes=1
#SBATCH --ntasks=12

#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=12
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

cd /scratch/alpine/noso3320/trem2/trem2_unidock
conda activate /projects/noso3320/hts_env
python $SCRIPT_DIR/dock.py
