#!/bin/bash
#SBATCH --job-name="test_CONUS"
#SBATCH -p cascade,sapphire
#SBATCH -t 12:00:00
#SBATCH --mem 4G
#SBATCH -o ../log/test_CONUS.out
#SBATCH -e ../log/test_CONUS.err

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate reactor-nu

python ../test_reactor_measurement.py -d ../configs/CONUS_detector.cfg -r ../configs/CONUS_reactor.cfg