#!/bin/bash
#SBATCH --partition=longq
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=20
#SBATCH --mem=60G
#SBATCH --exclude=compute28
#SBATCH --time=120:00:00
#SBATCH --output=log.barrier_%j.out
# ONE longq job on 8 nodes x 20 tasks (the longq allocation): lane 1 and 2 = 3 nodes each, lane 3 = the remaining 2 nodes. All lanes pull constrained-scan images
# (bar_*.in, listed in ads/list_bar_*.txt) from a shared queue, so Mg, Zn and later Ca are all handled by this single job.
# Before every claim the lists are rebuilt (13_barrier.py build, light, on the compute node), so Ca joins automatically once its
# dissociation rerun has converged. Images whose output has JOB DONE are skipped; resubmitting continues where it stopped.
module purge
module load apps/quantum-espresso/7.5/cpu
export OMP_NUM_THREADS=1
HOSTS=($(scontrol show hostnames "$SLURM_JOB_NODELIST"))
LOCK=ads/locks_$SLURM_JOB_ID; mkdir -p $LOCK
NL=(3 3 $(( ${#HOSTS[@]} - 6 )))      # nodes per lane
START=(0 3 6)
for L in 1 2 3; do : > hosts_lane$L.$SLURM_JOB_ID; done
for L in 1 2 3; do
  for ((i=0; i<${NL[$((L-1))]}; i++)); do echo "${HOSTS[$((START[$((L-1))]+i))]} slots=20" >> hosts_lane$L.$SLURM_JOB_ID; done
done
echo "=== barrier started: $(date) | nodes: ${HOSTS[@]} | per lane: ${NL[@]}"
refresh() { ( module load apps/anaconda3/2024.06 >/dev/null 2>&1; python3 13_barrier.py build > ads/barrier_build.log 2>&1 ); }
lane() {
  L=$1
  while true; do
    refresh
    NAME=""
    for LS in ads/list_bar_*.txt; do
      for N in $(cat $LS 2>/dev/null); do
        grep -q "JOB DONE" qe_inputs/$N.out 2>/dev/null && continue
        mkdir $LOCK/$N 2>/dev/null && { NAME=$N; break 2; }
      done
    done
    [ -z "$NAME" ] && break
    echo "lane$L: $NAME  $(date)"
    T0=$SECONDS
    mpirun -np $((20*${NL[$((L-1))]})) --hostfile hosts_lane$L.$SLURM_JOB_ID --map-by ppr:20:node pw.x -nk 1 -in qe_inputs/$NAME.in > qe_inputs/$NAME.out 2>&1
    echo "lane$L: $NAME -> $(grep -E 'bfgs converged|convergence NOT|Error' qe_inputs/$NAME.out | tail -2 | tr '\n' ' ') $(date)"
    if ! grep -q "JOB DONE" qe_inputs/$NAME.out && [ $((SECONDS-T0)) -lt 120 ]; then
      echo "lane$L: $NAME stopped after $((SECONDS-T0)) s without JOB DONE - aborting this lane. Last lines of the output:"
      tail -8 qe_inputs/$NAME.out; break
    fi
  done
  echo "lane$L: nothing left to claim"
}
lane 1 & sleep 20; lane 2 & sleep 20; [ ${NL[2]} -ge 1 ] && lane 3 &
wait
rm -f hosts_lane*.$SLURM_JOB_ID
echo "=== finished: $(date) ==="
