#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=log.prep_cref_%j.out
# Builds the site-matched clean-slab reference inputs from the existing ads_*.in files (compute node).
module load apps/anaconda3/2024.06
export CELLS="${CELLS:-crystal_ideal matrix_xp}"
python3 11_adsorption.py cleanref || exit 1
echo; ls qe_inputs | grep -c "^cref_"; wc -l ads/list_cref_*.txt
