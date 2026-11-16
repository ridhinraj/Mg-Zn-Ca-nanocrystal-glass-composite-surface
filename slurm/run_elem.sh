#!/bin/bash
# Element-stratified H2O sites on 5 new cells + water dissociation on 7 finished sites. Only sbatch calls: 6 jobs.
#   bash run_elem.sh        DRY=1 bash run_elem.sh
CELLS="${CELLS:-matrix_xm matrix_yp iface_S_D28 iface_N_D26 iface_S_D26}"; export CELLS
set -- $CELLS
sub() { if [ -n "$DRY" ]; then echo "DRY: sbatch $*" >&2; echo 0000; else sbatch --parsable "$@"; fi; }
P=$(sub prep_elem.sh); echo "prep: $P"
D="--dependency=afterok:$P --exclude=compute28 --time=80:00:00 --mem=60G"
pat()  { echo "ads_${1}_*_H2O"; }
lst()  { echo "ads/list_cref_${1}_a.txt ads/list_cref_${1}_b.txt"; }
J1=$(sub $D --export=ALL,PATS="$(pat $1)",LISTS="$(lst $1)" submit_ads.sh);                                       echo "$1: $J1"
J2=$(sub $D --export=ALL,PATS="$(pat $2)",LISTS="$(lst $2)" submit_ads.sh);                                       echo "$2: $J2"
J3=$(sub $D --export=ALL,PATS="$(pat $3) $(pat $4)",LISTS="$(lst $3) $(lst $4)" submit_ads.sh);                   echo "$3 + $4: $J3"
J4=$(sub $D --export=ALL,PATS="$(pat $5)",LISTS="$(lst $5) ads/list_diss_1.txt" submit_ads.sh);                   echo "$5 + dissociation 1: $J4"
J5=$(sub $D --export=ALL,LISTS="ads/list_diss_2.txt ads/list_diss_3.txt" submit_ads.sh);                          echo "dissociation 2 + 3: $J5"
echo "when all finished:  sbatch submit_ads_analyze.sh"
