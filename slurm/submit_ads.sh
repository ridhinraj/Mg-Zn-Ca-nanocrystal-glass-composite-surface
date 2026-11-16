#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --nodes=3
#SBATCH --ntasks=20
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#SBATCH --output=log.ads_%j.out
#SBATCH --error=log.ads_%j.err
# Runs every structure named in the list files (one per line), one after another:
#   sbatch --export=ALL,LISTS="ads/list_crystal.txt ads/list_hi.txt" submit_ads.sh
# Or name structures directly:  --export=ALL,NAMES="ads_matrix_xp_s4_H2O ads_crystal_ideal_atop_H2O_hi"
# Or by pattern (expanded when the job runs):  --export=ALL,PATS="ads_iface_S_D27_*_H2O cref_iface_S_D27_1"
# A structure whose output already has JOB DONE is skipped, so a resubmission continues where it stopped.
module purge
module load apps/quantum-espresso/7.5/cpu
export OMP_NUM_THREADS=1
echo "=== ads batch started: $(date) | $LISTS ==="
ALL=""
for L in $LISTS; do ALL="$ALL $(cat $L)"; done
ALL="$ALL $NAMES"
for P in $PATS; do for F in qe_inputs/$P.in; do [ -f "$F" ] && ALL="$ALL $(basename $F .in)"; done; done
for NAME in $ALL; do
  {
    OUT=qe_inputs/${NAME}.out
    if grep -q "JOB DONE" "$OUT" 2>/dev/null; then echo "skip $NAME (done)"; continue; fi
    [ -f qe_inputs/${NAME}.in ] || { echo "ERROR: missing qe_inputs/${NAME}.in"; continue; }
    echo "--- $NAME : $(date)"
    mpirun -np $SLURM_NTASKS pw.x -nk 1 -in qe_inputs/${NAME}.in > $OUT 2>&1
    grep -E "^!|bfgs converged|Final energy|JOB DONE|convergence NOT|Error" $OUT | tail -4
  }
done
echo "=== finished: $(date) ==="
