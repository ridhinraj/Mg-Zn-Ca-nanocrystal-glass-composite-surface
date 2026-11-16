#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=log.prep_barrier_refine_%j.out
# Extra images (lambda 0.22/0.26/0.34/0.38 for Mg, Zn; 0.38/0.42/0.50/0.54 for Ca) built from the relaxed neighbouring images.
module load apps/anaconda3/2024.06
python3 13_barrier.py refine
