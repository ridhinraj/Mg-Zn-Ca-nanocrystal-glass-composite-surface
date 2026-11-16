#!/bin/bash
# Pilot DFT set at Gamma: generate inputs, submit each QE job, and chain its projwfc job to start
# only if the QE job succeeded (--dependency=afterok).
#
#   bash run_pilot.sh          # submit for real
#   DRY=1 bash run_pilot.sh    # just print what would be submitted
#
# Run from W3/ (needs slab_<name>.xyz for each patch, the pseudopotentials, and the submit_*.sh files).

SMALL="crystal_ideal iface_N_D26 iface_S_D26"      # 360-460 electrons: minutes
BIG="matrix_xp matrix_xm"                          # 1300-1440 electrons, ~800 bands: hours

run() { if [ -n "$DRY" ]; then echo "DRY: $*"; else "$@"; fi; }
sub() { if [ -n "$DRY" ]; then echo "DRY: sbatch $*" >&2; echo 0000; else sbatch --parsable "$@"; fi; }

# 1. inputs (nbnd = occupied + 110, reaching ~E_F + 2 eV; Gamma only)
FILES=""; for n in $SMALL $BIG; do
    [ -f slab_$n.xyz ] || { echo "ERROR: slab_$n.xyz missing (run 07_build_periodic_patch.py export ...)"; exit 1; }
    FILES="$FILES slab_$n.xyz"
done
run python3 05_make_qe_inputs.py $FILES

# 2. small patches: default header of submit_qe.sh (mediumq, 3 nodes, 20 tasks, 40G per node)
for n in $SMALL; do
    jq=$(sub --export=ALL,NAME=$n --time=02:00:00 submit_qe.sh)
    jp=$(sub --dependency=afterok:$jq --export=ALL,NAME=$n submit_projwfc.sh)
    echo "$n : QE job $jq  ->  projwfc job $jp"
done

# 3. matrix patches: long limit; projwfc gets more time and memory (about 950 atomic wavefunctions)
for n in $BIG; do
    jq=$(sub --export=ALL,NAME=$n --time=12:00:00 submit_qe.sh)
    jp=$(sub --dependency=afterok:$jq --partition=mediumq --time=04:00:00 --mem=60G --export=ALL,NAME=$n submit_projwfc.sh)
    echo "$n : QE job $jq  ->  projwfc job $jp"
done

cat << 'MSG'

When the jobs have finished:
   for n in crystal_ideal iface_N_D26 iface_S_D26 matrix_xp matrix_xm; do python3 06_analyze_pdos.py $n; done
   sacct -u $USER -S today --format=JobID,JobName%20,State,Elapsed,ReqMem,MaxRSS | grep -v batch
and send me the printed analysis plus the QE summary lines (Estimated RAM, iterations, wall time).
MSG
