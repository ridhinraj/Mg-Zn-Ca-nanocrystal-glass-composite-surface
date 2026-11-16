#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --nodes=3
#SBATCH --ntasks=20
#SBATCH --mem=40G
#SBATCH --time=02:00:00
#SBATCH --output=log.batch_%j.out
#SBATCH --error=log.batch_%j.err
# ONE job per group of patches: for each NAME run pw.x, then projwfc.x, then the analysis,
# all inside this allocation. (The cluster limits how many jobs one user can have queued, so
# QE + projwfc + analysis are not separate jobs.)  Choose patches at submit time:
#   sbatch --export=ALL,NAMES="iface_N_D27 iface_S_D27" submit_batch.sh
# For matrix patches override:  --time=06:00:00 --mem=60G
NK="${NK:-1}"
module load apps/anaconda3/2024.06   # numpy for the analysis (compute-node /usr/bin/python3 has none)
PY="${PYBIN:-$(which python3)}"
echo "python: $PY"
$PY -c "import numpy" || { echo "ERROR: no numpy in $PY"; exit 1; }
module purge
module load apps/quantum-espresso/7.5/cpu
export OMP_NUM_THREADS=1
echo "=== batch started: $(date) | $NAMES ==="
for NAME in $NAMES; do
  echo; echo "################ $NAME : $(date)"
  INP=qe_inputs/${NAME}.in; OUT=qe_inputs/${NAME}.out
  [ -f "$INP" ] || { echo "ERROR: missing $INP"; continue; }
  mpirun -np $SLURM_NTASKS pw.x -nk $NK -in "$INP" > "$OUT" 2>&1
  grep -E "number of Kohn-Sham|convergence has been achieved|Fermi energy|JOB DONE" "$OUT"
  grep -q "JOB DONE" "$OUT" || { echo "ERROR: pw.x did not finish for $NAME - skipping"; continue; }

  EF=$(grep "the Fermi energy is" "$OUT" | tail -1 | awk '{print $5}')
  EMIN=$(awk -v e=$EF 'BEGIN{printf "%.2f", e-12.0}'); EMAX=$(awk -v e=$EF 'BEGIN{printf "%.2f", e+4.0}')
  mkdir -p pdos/${NAME}
  cat > qe_inputs/${NAME}.projwfc.in << EOT
&projwfc
  prefix  = '${NAME}'
  outdir  = './tmp_${NAME}'
  filpdos = 'pdos/${NAME}/${NAME}'
  Emin    = ${EMIN}
  Emax    = ${EMAX}
  DeltaE  =  0.02
  ngauss  =  0
  degauss =  0.01
/
EOT
  mpirun -np 4 projwfc.x -nk 1 -in qe_inputs/${NAME}.projwfc.in > qe_inputs/${NAME}.projwfc.out 2>&1
  grep -q "JOB DONE" qe_inputs/${NAME}.projwfc.out || { echo "ERROR: projwfc failed for $NAME"; continue; }
  echo "--- analysis $NAME ---"
  $PY 06_analyze_pdos.py $NAME
done
echo "=== batch finished: $(date) ==="
