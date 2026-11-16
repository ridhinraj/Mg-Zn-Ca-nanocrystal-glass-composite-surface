#!/bin/bash
# Barrier scans with the proven submit_ads.sh (nodes x 20 tasks TOTAL, plain mpirun - no hostfile). Only sbatch calls.
# Lists ads/list_bar_{Mg,Zn,Ca}.txt already exist from prep_barrier.sh (re-run  sbatch prep_barrier.sh  if they are missing).
#   Mg, Zn : mediumq, 3 nodes, 20 tasks          Ca : longq, 8 nodes, 20 tasks (2-3 tasks per node, easy to schedule)
#   bash run_barrier.sh       DRY=1 bash run_barrier.sh
sub() { if [ -n "$DRY" ]; then echo "DRY: sbatch $*" >&2; echo 0000; else sbatch --parsable "$@"; fi; }
C="--exclude=compute28 --mem=60G --time=80:00:00"
for EL in Mg Zn; do
  J=$(sub --partition=mediumq --nodes=3 $C --export=ALL,LISTS="ads/list_bar_$EL.txt" submit_ads.sh); echo "$EL (mediumq): $J"
done
J=$(sub --partition=longq --nodes=8 $C --export=ALL,LISTS="ads/list_bar_Ca.txt" submit_ads.sh); echo "Ca (longq): $J"
echo "When finished:  sbatch submit_barrier_analyze.sh"
