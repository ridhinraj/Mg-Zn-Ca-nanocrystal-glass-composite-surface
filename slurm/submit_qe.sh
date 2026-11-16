#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --nodes=3
#SBATCH --ntasks=20
#SBATCH --mem=40G
#SBATCH --time=02:00:00
#SBATCH --output=log.qe_%j.out
#SBATCH --error=log.qe_%j.err

# LAYOUT / MEMORY NOTES
#  * --mem is PER NODE and needs a unit (40G, not 40 = 40 MB).
#  * QE's printed estimate is ~63 GB TOTAL (at the old 320 Ry cutoff; less at 200), split
#    over the nodes, so 3 nodes -> ~21 GB/node (40G gives 2x margin), 8 nodes -> ~8 GB/node.
#  * Layouts you can try (pick whichever starts soonest; speed scales roughly with ranks):
#      mediumq, 3 nodes, 20 tasks            : as set above (~22 s/iteration measured, 165-atom Mg patch)
#      mediumq, 3 nodes, 10 tasks per node   : --nodes=3 --ntasks-per-node=10 (30 tasks), --mem=40G
#      longq,   8 nodes, 40 tasks            : --partition=longq --nodes=8 --ntasks=40 --mem=20G
#  * See which would start first WITHOUT submitting:
#      sbatch --test-only --partition=mediumq --nodes=3 --ntasks=20 --mem=40G --time=02:00:00 submit_qe.sh
#  * Override at submit time with e.g.  sbatch --nodes=8 --ntasks=40 --partition=longq --mem=20G submit_qe.sh

# One patch per job.  Default = pure-Mg crystal core (338 e-, cheapest: a real
# end-to-end test).  Choose another patch at submit time with e.g.
#   sbatch --export=ALL,NAME=Mg60Zn35Ca5_interphase submit_qe.sh
# For the big patches use a longer partition and more cores, e.g.
#   sbatch --partition=shortq --nodes=2 --ntasks=80 --time=24:00:00 \
#          --export=ALL,NAME=Mg60Zn35Ca5_glass_core submit_qe.sh

NAME="${NAME:-Mg60Zn35Ca5_crystal_core}"
NK="${NK:-1}"        # k-point pools. Use NK=1 for Gamma-only; NK=2 for the 2x2x1 / 3x3x1 tests (memory grows ~x NK)
INP="qe_inputs/${NAME}.in"
OUT="qe_inputs/${NAME}.out"

echo "=== QE job started: $(date) | node: $(hostname) | patch: $NAME | pools: $NK ==="

# Clear the login-node default MPI, then load QE (brings GCC 10.5 + OpenMPI 4.1.8 + MKL).
module purge
module load apps/quantum-espresso/7.5/cpu
module list 2>&1
echo "pw.x:   $(which pw.x)"
echo "mpirun: $(which mpirun)"

export OMP_NUM_THREADS=1

if [ ! -f "$INP" ]; then
    echo "ERROR: missing $INP  (run: python3 05_make_qe_inputs.py)"; exit 1
fi
PSEUDO=$(grep -m1 pseudo_dir "$INP" | cut -d"'" -f2)
if [ ! -d "$PSEUDO" ]; then
    echo "ERROR: pseudo_dir $PSEUDO not found"; exit 1
fi

mpirun -np $SLURM_NTASKS pw.x -nk $NK -in "$INP" > "$OUT" 2>&1

echo "=== QE job finished: $(date) ==="
echo "--- after the job ends, run:  seff $SLURM_JOB_ID   and   sacct -j $SLURM_JOB_ID --format=JobID,State,ReqMem,MaxRSS,NodeList ---"
echo "--- summary from $OUT ---"
grep -i -E "writing|running on|MPI processes|Estimated|number of electrons|number of Kohn-Sham|convergence has been achieved|^!|Fermi energy|Total force|JOB DONE" "$OUT"
if ! grep -q "JOB DONE" "$OUT"; then
    echo "WARNING: no 'JOB DONE' - check the end of $OUT"
    exit 1          # non-zero exit so a chained projwfc job (--dependency=afterok) does not start
fi
