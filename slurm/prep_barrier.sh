#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=log.prep_barrier_%j.out
# Builds the constrained-scan inputs and lists ads/list_bar_<El>.txt for the elements that have a converged dissociation pair.
module load apps/anaconda3/2024.06
python3 13_barrier.py build
