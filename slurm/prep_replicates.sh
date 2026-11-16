#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=log.prep_%j.out
#SBATCH --error=log.prep_%j.err
# Compute-node job: build 6 replicate patches, MEAM-relax them, export xyz, write QE inputs.
# Nothing here runs on the login node. Submit with:  bash run_replicates.sh
PARENT=AN_composite_relaxed_300K.data
MATRIX="matrix_yp matrix_ym"
IFACE="iface_N_D27 iface_S_D27 iface_N_D28 iface_S_D28"
ALL="$MATRIX $IFACE"
set -e
module load apps/anaconda3/2024.06   # compute nodes' /usr/bin/python3 has no numpy
PY="${PYBIN:-$(which python3)}"      # absolute path, survives later module loads
echo "python: $PY"
$PY -c "import numpy" || { echo "ERROR: $PY has no numpy on $(hostname)"; exit 1; }
echo "=== prep started: $(date) ==="
# python (numpy only) - same environment as the earlier jobs; edit the module line if yours differs
$PY 07_build_periodic_patch.py build --kind matrix --offset 0 60  --name matrix_yp --parent $PARENT
$PY 07_build_periodic_patch.py build --kind matrix --offset 0 -60 --name matrix_ym --parent $PARENT
for D in 27 28; do
  $PY 07_build_periodic_patch.py build --kind interphase --D $D --pole north --name iface_N_D$D --parent $PARENT
  $PY 07_build_periodic_patch.py build --kind interphase --D $D --pole south --name iface_S_D$D --parent $PARENT
done
set +e
module load apps/lammps/intel/parallel_studio_xe_2020.4.912/22Jul2025/Optional/1
LAMMPS_BIN=$HOME/lammps_meam/bin/lmp
MPI_BIN=/opt/apps/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpirun
for NAME in $ALL; do
  RB=0; case "$NAME" in matrix*) RB=1;; esac
  echo "--- relax $NAME (relax_box=$RB) ---"
  $MPI_BIN -np $SLURM_NTASKS $LAMMPS_BIN -in 08_relax_patch.in \
     -var patch_data patches/${NAME}.data -var out patches/${NAME} -var relax_box $RB \
     -log patches/log.relax_${NAME} > /dev/null 2>&1
  grep -E "^FINAL" patches/log.relax_${NAME} || { echo "ERROR: relax of $NAME failed"; exit 1; }
done
for n in $ALL; do $PY 07_build_periodic_patch.py export patches/${n}_relaxed.data || exit 1; done
$PY 05_make_qe_inputs.py $(for n in $ALL; do echo slab_$n.xyz; done) || exit 1
echo "=== prep finished: $(date) ==="
