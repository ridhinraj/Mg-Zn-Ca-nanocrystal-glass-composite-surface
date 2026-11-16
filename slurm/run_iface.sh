#!/bin/bash
# Interphase adsorption (H2O and OH) with site-matched clean references. Only sbatch calls: 6 jobs.
#   bash run_iface.sh        DRY=1 bash run_iface.sh      CELLS="iface_S_D27 iface_N_D28" (default; one cell per pole)
CELLS="${CELLS:-iface_S_D27 iface_N_D28}"; export CELLS
set -- $CELLS; C1=$1; C2=$2
sub() { if [ -n "$DRY" ]; then echo "DRY: sbatch $*" >&2; echo 0000; else sbatch --parsable "$@"; fi; }
P=$(sub prep_iface.sh); echo "prep: $P"
D="--dependency=afterok:$P --exclude=compute28 --time=60:00:00 --mem=60G"
J1=$(sub $D --export=ALL,PATS="ads_${C1}_*_H2O",LISTS="ads/list_cref_${C1}_a.txt" submit_ads.sh);  echo "$C1 H2O + references a: $J1"
J2=$(sub $D --export=ALL,PATS="ads_${C1}_*_OH"  submit_ads.sh);                                    echo "$C1 OH: $J2"
J3=$(sub $D --export=ALL,PATS="ads_${C2}_*_H2O",LISTS="ads/list_cref_${C2}_a.txt" submit_ads.sh);  echo "$C2 H2O + references a: $J3"
J4=$(sub $D --export=ALL,PATS="ads_${C2}_*_OH"  submit_ads.sh);                                    echo "$C2 OH: $J4"
J5=$(sub $D --export=ALL,LISTS="ads/list_cref_${C1}_b.txt ads/list_cref_${C2}_b.txt" submit_ads.sh); echo "references b (both cells): $J5"
echo "when all finished:  sbatch submit_ads_analyze.sh"
