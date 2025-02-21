#!/bin/bash
#SBATCH --job-name="RELICS_run1"
#SBATCH -p cascade,sapphire
#SBATCH -t 12:00:00
#SBATCH --mem 4G
#SBATCH -o ../log/RELICS_run1.out
#SBATCH -e ../log/RELICS_run1.err

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate reactor-nu

python ../run_reactor_measurement.py -d ../configs/RELICS_detector.cfg -r ../configs/RELICS_reactor.cfg -o RELICS_run1
