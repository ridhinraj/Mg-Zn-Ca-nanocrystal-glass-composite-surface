#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=log.barrier_analyze_%j.out
module load apps/anaconda3/2024.06
python3 13_barrier.py analyze
