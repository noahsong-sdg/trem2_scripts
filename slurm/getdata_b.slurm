#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=preemptable
#SBATCH --job-name=get_data_b
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=datab_%j.out
#SBATCH --error=datab_%j.err

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

python scripts/getdatab.py # c2
python scripts/getdatab2.py # c3
#python scripts/getdata3.py # c1
