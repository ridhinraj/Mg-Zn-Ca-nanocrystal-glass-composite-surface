#!/bin/bash
# Only sbatch calls: 5 jobs (prep, crystal, 3 matrix lists). Needs the finished adsorption run in place.
#   bash run_cref.sh        DRY=1 bash run_cref.sh
sub() { if [ -n "$DRY" ]; then echo "DRY: sbatch $*" >&2; echo 0000; else sbatch --parsable "$@"; fi; }
P=$(sub prep_cref.sh);                                                                              echo "prep: $P"
D="--dependency=afterok:$P --exclude=compute28"
C=$(sub $D --time=04:00:00 --export=ALL,LISTS="ads/list_cref_crystal.txt" submit_ads.sh);            echo "crystal clean refs: $C"
for k in 1 2 3; do
  M=$(sub $D --time=30:00:00 --mem=60G --export=ALL,LISTS="ads/list_cref_matrix_$k.txt" submit_ads.sh); echo "matrix clean refs $k: $M"
done
echo "when all finished:  sbatch submit_ads_analyze.sh"
