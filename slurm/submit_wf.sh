#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --nodes=3
#SBATCH --ntasks=20
#SBATCH --mem=40G
#SBATCH --time=03:00:00
#SBATCH --output=log.wf_%j.out
#SBATCH --error=log.wf_%j.err
# Work function: for each NAME build the dipole-corrected input, run pw.x, pp.x (potential), then analyse.
#   sbatch --export=ALL,NAMES="crystal_ideal iface_N_D26" submit_wf.sh                       # validation first
#   sbatch --time=06:00:00 --mem=60G --export=ALL,NAMES="matrix_xp matrix_xm matrix_yp matrix_ym" submit_wf.sh
module load apps/anaconda3/2024.06
PY=$(which python3)
module purge
module load apps/quantum-espresso/7.5/cpu
export OMP_NUM_THREADS=1
mkdir -p wf
echo "=== work-function batch started: $(date) | $NAMES ==="
$PY 09_make_wf_inputs.py $NAMES || exit 1
for NAME in $NAMES; do
  echo; echo "################ $NAME : $(date)"
  mpirun -np $SLURM_NTASKS pw.x -nk 1 -in qe_inputs/${NAME}_wf.in > qe_inputs/${NAME}_wf.out 2>&1
  grep -E "convergence has been achieved|Fermi energy|JOB DONE|Dipole|dipole" qe_inputs/${NAME}_wf.out
  grep -q "JOB DONE" qe_inputs/${NAME}_wf.out || { echo "ERROR: pw.x failed for $NAME"; continue; }
  cat > wf/${NAME}.pp.in << EOT
&inputpp
  prefix = '${NAME}_wf'
  outdir = './tmp_${NAME}_wf'
  filplot = 'wf/${NAME}.pot'
  plot_num = 11
/
&plot
  iflag = 3
  output_format = 6
  fileout = 'wf/${NAME}.cube'
/
EOT
  mpirun -np 1 pp.x -in wf/${NAME}.pp.in > wf/${NAME}.pp.out 2>&1
  [ -s wf/${NAME}.cube ] || { echo "ERROR: pp.x produced no cube for $NAME (see wf/${NAME}.pp.out)"; continue; }
  $PY 10_workfunction.py $NAME
done
echo "=== finished: $(date) ==="
