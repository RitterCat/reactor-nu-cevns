#!/bin/bash
#SBATCH --job-name="breeder_1"
#SBATCH -p sapphire
#SBATCH -t 12:00:00
#SBATCH --mem 4G
#SBATCH -o ../log/PWR_breeder_run1.out
#SBATCH -e ../log/PWR_breeder_run1.err

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate reactor-nu

python ../run_reactor_measurement.py -d ../configs/germanium_detector_no_exposure.cfg --exposure 1 -r ../configs/PWR_breeder.cfg -o PWR_breeder_run1
