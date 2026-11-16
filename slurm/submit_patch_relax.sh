#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=log.patch_relax_%j.out
#SBATCH --error=log.patch_relax_%j.err

# Relax one or more periodic patches (170 atoms each: minutes, not hours).
#   sbatch --export=ALL,NAMES="crystal_ideal iface_N_D26 iface_S_D26 matrix_xp matrix_xm" submit_patch_relax.sh
# Names starting with "matrix" also relax the lateral cell (amorphous); the rest keep it.
NAMES="${NAMES:-crystal_ideal iface_N_D26 iface_S_D26 matrix_xp matrix_xm}"

echo "=== patch relax started: $(date) | $NAMES ==="
module load apps/lammps/intel/parallel_studio_xe_2020.4.912/22Jul2025/Optional/1
LAMMPS_BIN=$HOME/lammps_meam/bin/lmp
MPI_BIN=/opt/apps/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpirun

for f in library.meam CaZnMg.meam 08_relax_patch.in; do
    [ -f "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done

for NAME in $NAMES; do
    [ -f patches/${NAME}.data ] || { echo "SKIP $NAME: patches/${NAME}.data not found"; continue; }
    RB=0; case "$NAME" in matrix*) RB=1;; esac
    echo "--- $NAME (relax_box=$RB) ---"
    $MPI_BIN -np $SLURM_NTASKS $LAMMPS_BIN -in 08_relax_patch.in \
        -var patch_data patches/${NAME}.data -var out patches/${NAME} -var relax_box $RB \
        -log patches/log.relax_${NAME} > /dev/null 2>&1
    grep -E "^FINAL" patches/log.relax_${NAME} || echo "WARNING: $NAME did not finish - see patches/log.relax_${NAME}"
done
echo "=== patch relax finished: $(date) ==="
