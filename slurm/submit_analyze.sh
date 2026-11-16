#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=log.analyze_%j.out
# Per-atom analysis on a compute node.
#   sbatch --export=ALL,NAMES="matrix_yp matrix_ym iface_N_D27 iface_S_D27 iface_N_D28 iface_S_D28" submit_analyze.sh
module load apps/anaconda3/2024.06   # compute-node /usr/bin/python3 has no numpy
for n in $NAMES; do echo "##### $n"; python3 06_analyze_pdos.py $n; done
