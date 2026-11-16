#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=log.prep_ads_%j.out
# Builds sites + relax inputs + gas references (compute node). Needs analysis_matrix_xp.csv and
# qe_inputs/{crystal_ideal,matrix_xp}_wf.in in the submit directory.
module load apps/anaconda3/2024.06
PY=$(which python3)
echo "python: $PY"
export CELLS="${CELLS:-crystal_ideal matrix_xp}"
$PY 11_adsorption.py build || exit 1
echo; ls qe_inputs | grep -c "^ads_" ; ls ads/
