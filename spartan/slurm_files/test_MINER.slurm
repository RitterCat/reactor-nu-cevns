#!/bin/bash
#SBATCH --job-name="test_MINER"
#SBATCH -p cascade,sapphire
#SBATCH -t 12:00:00
#SBATCH --mem 4G
#SBATCH -o ../log/test_MINER.out
#SBATCH -e ../log/test_MINER.err

module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"
conda activate reactor-nu

python ../test_reactor_measurement.py -d ../configs/MINER_detector.cfg -r ../configs/HFIR.cfg