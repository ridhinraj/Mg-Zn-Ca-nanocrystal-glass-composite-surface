#!/bin/bash
# Composition series for the merged corrosion paper: cut + relax 2 glass-matrix patches per composition,
# export them and build the QE inputs the surface model reads. NO DFT. Only sbatch calls.
#
#   1. edit SERIES below: one "name:path" per composition (paths = Study 2's quenched matrix .data files)
#   2. bash run_series.sh            DRY=1 bash run_series.sh    to see the jobs first
#   3. when it finishes:  sbatch submit_series_model.sh
#
# Each entry's name becomes the patch prefix, so keep them short and sortable, e.g.
#   c04_zn25  c09_zn32  c01_zn35  c10_zn38   (Zn series at x_Ca = 0.05)
#   c05_ca03  c01_ca05  c06_ca07             (Ca series at x_Zn = 0.35; c01 is shared, cut once)
SERIES="${SERIES:-\
c04_zn25:/CHANGE/ME/comp4.data \
c09_zn32:/CHANGE/ME/comp9.data \
c01_zn35:/CHANGE/ME/comp1.data \
c10_zn38:/CHANGE/ME/comp10.data \
c05_ca03:/CHANGE/ME/comp5.data \
c06_ca07:/CHANGE/ME/comp6.data}"
OFFSETS="${OFFSETS:-60,0 0,60}"      # two independent, non-overlapping glass patches per composition
sub() { if [ -n "$DRY" ]; then echo "DRY: sbatch $*" >&2; echo 0000; else sbatch --parsable "$@"; fi; }
for e in $SERIES; do
  case "$e" in *:*) ;; *) echo "SKIP '$e' (expected name:path)"; continue;; esac
done
export SERIES OFFSETS
J=$(sub prep_series.sh); echo "cut + relax + export + QE inputs: $J"
echo "when it finishes:  sbatch submit_series_model.sh"
