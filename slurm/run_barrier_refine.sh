#!/bin/bash
# Prep (testq) then one job per element on the extra images. Mg, Zn: mediumq 3 nodes x 20 tasks; Ca: longq 8 nodes x 20 tasks.
#   bash run_barrier_refine.sh        DRY=1 bash run_barrier_refine.sh
sub() { if [ -n "$DRY" ]; then echo "DRY: sbatch $*" >&2; echo 0000; else sbatch --parsable "$@"; fi; }
P=$(sub prep_barrier_refine.sh); echo "prep: $P"
C="--dependency=afterok:$P --exclude=compute28 --mem=60G --time=80:00:00"
for EL in Mg Zn; do
  J=$(sub --partition=mediumq --nodes=3 $C --export=ALL,LISTS="ads/list_bar_ref_$EL.txt" submit_ads.sh); echo "$EL (mediumq): $J"
done
J=$(sub --partition=longq --nodes=8 $C --export=ALL,LISTS="ads/list_bar_ref_Ca.txt" submit_ads.sh); echo "Ca (longq): $J"
echo "When finished:  sbatch submit_barrier_analyze.sh"
