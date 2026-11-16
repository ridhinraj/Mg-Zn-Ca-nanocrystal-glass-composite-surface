#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=log.pband_%j.out
# PYTHONNOUSERSITE=1 hides ~/.local (numpy 2.2.6) so anaconda's own numpy 1.x is used and its
# matplotlib imports. Everything here works on numpy 1.x.
module load apps/anaconda3/2024.06
export PYTHONNOUSERSITE=1
python3 -c "import numpy,matplotlib;print('numpy',numpy.__version__,'| matplotlib',matplotlib.__version__)" || {
  echo "ERROR: matplotlib still unavailable with PYTHONNOUSERSITE=1"; exit 1; }
python3 23_pband_region.py
