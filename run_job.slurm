#!/bin/bash

#SBATCH --job-name=dummy                              # Job name
#SBATCH --ntasks=40                                   # Number of CPUs
#SBATCH --time=168:00:00                              # Time limit hrs:min:sec
#SBATCH --output=log_files/dummy_slurm.log            # Standard output and error log
#SBATCH --exclude=gpu01                 # Don't hog the GPU node


# memory
ulimit -s unlimited

# location
cd /home/gsnearing/projects/multi-inputs/SACSMA-SNOW17 
pwd

# load conda environment
source /home/gsnearing/.bashrc
conda activate sacsma

# run code
python batch_run_experiment.py --config_file=${1} --max_model_runs=${2} --algorithm=${3} --out_dir=${4}
