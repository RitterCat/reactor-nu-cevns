#!/bin/bash
#SBATCH --job-name="Ge_run2"
#SBATCH -p cascade,sapphire
#SBATCH -t 12:00:00
#SBATCH --mem 4G
#SBATCH -o ../log/germanium_run2.out
#SBATCH -e ../log/germanium_run2.err

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate reactor-nu

python ../run_reactor_measurement.py -d ../configs/germanium_detector.cfg -r ../configs/PWR_reactor.cfg -o germanium_run2 --power_prior
