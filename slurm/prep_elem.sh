#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=log.prep_elem_%j.out
# (1) element-stratified sites + H2O inputs + site-matched clean references for new cells
# (2) water-dissociation inputs (H2O* -> OH* + H*) for already finished H2O sites
# Needs qe_inputs/<cell>_wf.in and analysis_<cell>.csv for every cell in CELLS, and the finished ads_*_H2O outputs.
module load apps/anaconda3/2024.06
export CELLS="${CELLS:-matrix_xm matrix_yp iface_S_D28 iface_N_D26 iface_S_D26}"
export SITE_MODE=element MOLS=H2O N_SITES=5
python3 11_adsorption.py build    || exit 1
python3 11_adsorption.py cleanref || exit 1
python3 11_adsorption.py dissoc   || exit 1
echo; for c in $CELLS; do echo "$c: $(ls qe_inputs | grep -c "^ads_${c}_") H2O inputs, $(cat ads/list_cref_${c}_a.txt ads/list_cref_${c}_b.txt | grep -c .) clean references"; done
echo; cut -d, -f1-4 ads/sites_*.csv | grep -v "^site" | awk -F, '{print $3}' | sort | uniq -c
