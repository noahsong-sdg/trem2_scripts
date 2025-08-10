#SBATCH --time=1-00:00:00
#SBATCH --qos=preemptable
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --job-name=Ig
#SBATCH --output=mcdock.%j.out
#SBATCH --error=mcdock_%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=noso3320@colorado.edu

module purge
module load gcc
module load cuda/11.8
module load cmake
module load anaconda

eval "$(conda shell.bash hook)"
SCRIPT_DIR="/scratch/alpine/noso3320/trem2/unidock/scripts"
cd /scratch/alpine/noso3320/trem2/unidock
conda activate dock_env
python $SCRIPT_DIR/mcdock_b.py
