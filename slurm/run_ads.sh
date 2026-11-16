#!/bin/bash
# Only sbatch calls: 5 jobs (the cluster allows about 6 queued per user).
#   bash run_ads.sh        DRY=1 bash run_ads.sh
sub() { if [ -n "$DRY" ]; then echo "DRY: sbatch $*" >&2; echo 0000; else sbatch --parsable "$@"; fi; }
P=$(sub prep_ads.sh);                                                                                         echo "prep: $P"
D="--dependency=afterok:$P"
G=$(sub $D --partition=testq --nodes=1 --ntasks=8 --mem=16G --time=01:00:00 --export=ALL,LISTS="ads/list_gas.txt" submit_ads.sh);  echo "gas references: $G"
C=$(sub $D --time=08:00:00 --export=ALL,LISTS="ads/list_crystal.txt ads/list_hi.txt" submit_ads.sh);          echo "crystal_ideal + cutoff test: $C"
M1=$(sub $D --time=24:00:00 --mem=60G --export=ALL,LISTS="ads/list_matrix_H2O.txt" submit_ads.sh);            echo "matrix H2O: $M1"
M2=$(sub $D --time=24:00:00 --mem=60G --export=ALL,LISTS="ads/list_matrix_OH.txt" submit_ads.sh);             echo "matrix OH: $M2"
echo "when all finished:  sbatch submit_ads_analyze.sh"
