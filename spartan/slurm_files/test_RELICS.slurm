#!/bin/bash
#SBATCH --job-name="test_RELICS"
#SBATCH -p cascade,sapphire
#SBATCH -t 12:00:00
#SBATCH --mem 4G
#SBATCH -o ../log/test_RELICS.out
#SBATCH -e ../log/test_RELICS.err

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate reactor-nu

python ../test_reactor_measurement.py -d ../configs/RELICS_detector.cfg -r ../configs/RELICS_reactor.cfg