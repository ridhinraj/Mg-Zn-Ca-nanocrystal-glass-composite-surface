#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=log.prep_barrier_cont_%j.out
# Restart scan images that hit the 80-step BFGS limit (nstep=200 from their last geometry).
module load apps/anaconda3/2024.06
python3 13_barrier.py cont
