#!/bin/bash
# Only sbatch calls. 4 jobs total (cluster submit limit): prep, interface batch, 2 matrix batches.
#   bash run_replicates.sh            DRY=1 bash run_replicates.sh
# If the prep job is already queued/finished, pass its id:  PREP=117026 bash run_replicates.sh
sub() { if [ -n "$DRY" ]; then echo "DRY: sbatch $*" >&2; echo 0000; else sbatch --parsable "$@"; fi; }
if [ -z "$PREP" ]; then PREP=$(sub prep_replicates.sh); fi
echo "prep: $PREP"
D="--dependency=afterok:$PREP"
j1=$(sub $D --time=03:00:00 --export=ALL,NAMES="iface_N_D27 iface_S_D27 iface_N_D28 iface_S_D28" submit_batch.sh)
echo "interface batch: $j1 -> log.batch_$j1.out"
for n in matrix_yp matrix_ym; do
  j=$(sub $D --time=06:00:00 --mem=60G --export=ALL,NAMES="$n" submit_batch.sh)
  echo "$n : $j -> log.batch_$j.out"
done
