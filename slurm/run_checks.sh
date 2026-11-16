#!/bin/bash
# Tier-1 checks: one mediumq job for the cheap ones, one for the 60/480 site (the expensive one).
#   bash run_checks.sh        DRY=1 bash run_checks.sh
sub() { if [ -n "$DRY" ]; then echo "DRY: sbatch $*" >&2; echo 0000; else sbatch --parsable "$@"; fi; }
P=$(sub prep_checks.sh); echo "prep: $P"
C="--dependency=afterok:$P --exclude=compute28 --mem=60G --time=40:00:00 --nodes=3"
J1=$(sub --partition=mediumq $C --export=ALL,LISTS="ads/list_chk_kpt.txt ads/list_chk_wffix.txt ads/list_chk_d3.txt" submit_ads.sh)
echo "k-points + work function + dispersion: $J1"
J2=$(sub --partition=mediumq $C --export=ALL,LISTS="ads/list_chk_cut.txt" submit_ads.sh)
echo "60/480 cutoff site: $J2"
echo "when both finish:  sbatch submit_checks_analyze.sh"
