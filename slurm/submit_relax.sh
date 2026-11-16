#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --nodes=3
#SBATCH --ntasks=20
#SBATCH --time=7-00:00:00
#SBATCH --output=log.relax_surface_%j.out
#SBATCH --error=log.relax_surface_%j.err

# NOTE ON RESOURCES: unlike step 01 (geometry-only), this step runs real
# MD - 20,000 steps NVT anneal + 10,000 steps quench + minimization on
# ~194,619 atoms with pair_style meam. 6 hours is a generous first guess;
# check the log's reported loop time per step after it finishes and
# shrink the request for the next 2 Tier-1 compositions if it finished
# well under that.

echo "=== Surface Relaxation Job Started: $(date) ==="

module load apps/lammps/intel/parallel_studio_xe_2020.4.912/22Jul2025/Optional/1

LAMMPS_BIN=$HOME/lammps_meam/bin/lmp
MPI_BIN=/opt/apps/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpirun

# ---------------------------------------------------------------- INPUTS
# SLAB defaults to step 01's standard output; override with
# --export=SLAB=... when relaxing a different composition's cut slab.
SLAB="${SLAB:-dualphase_equatorial_slab_raw.data}"

echo "=== Node: $(hostname) | slab=$SLAB ==="

# 1. Check if the cut slab exists
if [ ! -f "$SLAB" ]; then
    echo "ERROR: Missing input file $SLAB - did step 01 (cut) complete?"
    exit 1
fi

# 2. Check potential files are present alongside it
for POTFILE in library.meam CaZnMg.meam; do
    if [ ! -f "$POTFILE" ]; then
        echo "ERROR: Missing potential file $POTFILE"
        exit 1
    fi
done

mkdir -p runs/corrosion_stage1

# 3. Run the surface relaxation
echo "--- Relaxing surface from $SLAB ---"
$MPI_BIN -np $SLURM_NTASKS $LAMMPS_BIN -in 02_relax_surface.in \
    -var slab_file $SLAB \
    -log runs/corrosion_stage1/log.02_relax_surface

echo "=== Surface Relaxation Job Finished: $(date) ==="
ls -la equatorial_slab_relaxed.data equatorial_slab_relaxed.dump 2>&1 || \
    echo "WARNING: expected output files not found - check runs/corrosion_stage1/log.02_relax_surface for errors"
