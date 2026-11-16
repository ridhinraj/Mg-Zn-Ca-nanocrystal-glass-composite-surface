#!/bin/bash
# OH at every finished H2O site: 1 prep job on testq + 3 DFT jobs on mediumq (leaves 1 mediumq slot free).
#   bash run_oh.sh            DRY=1 bash run_oh.sh         MAX_SITES=12 bash run_oh.sh
export NJOBS=3
sub() { if [ -n "$DRY" ]; then echo "DRY: sbatch $*" >&2; echo 0000; else sbatch --parsable "$@"; fi; }
P=$(sub --export=ALL,NJOBS=3,MAX_SITES="${MAX_SITES:-0}",CELLS="$CELLS" prep_oh.sh); echo "prep: $P"
D="--dependency=afterok:$P --exclude=compute28 --time=80:00:00 --mem=60G"
for k in 1 2 3; do
  J=$(sub $D --export=ALL,LISTS="ads/list_oh_${k}.txt" submit_ads.sh); echo "  OH batch $k: $J"
done
echo
echo "check the plan first:   cat log.prep_oh_${P}.out"
echo "when all finished:      sbatch submit_oh_report.sh"
