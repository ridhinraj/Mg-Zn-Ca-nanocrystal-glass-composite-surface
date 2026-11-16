#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=log.prep_iface_%j.out
# Builds sites, adsorption inputs and site-matched clean-reference inputs for the interphase cells (compute node).
# Needs qe_inputs/<cell>_wf.in and analysis_<cell>.csv for every cell in CELLS.
module load apps/anaconda3/2024.06
export CELLS="${CELLS:-iface_S_D27 iface_N_D28}"
python3 11_adsorption.py build    || exit 1
python3 11_adsorption.py cleanref || exit 1
echo; ls ads/; for c in $CELLS; do echo "$c: $(ls qe_inputs | grep -c "^ads_${c}_") adsorption inputs, $(cat ads/list_cref_${c}_a.txt ads/list_cref_${c}_b.txt | grep -c .) clean references"; done
