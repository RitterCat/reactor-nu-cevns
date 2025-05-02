#!/bin/bash
#SBATCH --job-name="brd_16"
#SBATCH -p sapphire
#SBATCH -t 12:00:00
#SBATCH --mem 4G
#SBATCH -o ../log/PWR_breeder_run16.out
#SBATCH -e ../log/PWR_breeder_run16.err

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate reactor-nu

python ../run_reactor_measurement.py -d ../configs/germanium_detector_no_exposure.cfg --exposure 0.8 -r ../configs/PWR_breeder.cfg --sampling_efficiency 0.3 -o PWR_breeder_run16