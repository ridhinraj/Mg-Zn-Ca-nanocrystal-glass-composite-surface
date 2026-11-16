#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=log.prep_oh_%j.out
# OH inputs at the existing H2O sites, plus any clean references they need. Additive: re-picks nothing,
# renumbers nothing, and must NOT be replaced by "11_adsorption.py build/cleanref" on these cells.
module load apps/anaconda3/2024.06
export NJOBS="${NJOBS:-3}"
python3 21_oh.py build || exit 1
