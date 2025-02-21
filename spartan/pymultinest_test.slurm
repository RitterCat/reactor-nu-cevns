#!/bin/bash
#SBATCH --job-name="pymultinest_test"
#SBATCH -p cascade,sapphire
#SBATCH -t 120
#SBATCH --mem 4G
#SBATCH -o ./log/test.out
#SBATCH -e ./log/test.err

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate reactor-nu

python pymultinest_test_script.py