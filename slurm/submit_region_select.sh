#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --output=log.region_select_%j.out
#SBATCH --error=log.region_select_%j.err

# Serial Python job (numpy only - no scipy, no OVITO, no MPI, no internet needed).
# Reads the relaxed slab dump (~194,619 atoms) and writes the three patch
# .xyz files. Should finish in about a minute.

echo "=== Region Selection Job Started: $(date) ==="
echo "=== Node: $(hostname) ==="

# TODO: if numpy/scipy only come from a module on your cluster, load it here,
# e.g.  module load apps/python/3.x   (check with: module avail python)
# module load apps/python/3.x

DUMP="${DUMP:-equatorial_slab_relaxed.dump}"

# 1. Check the relaxed slab exists
if [ ! -f "$DUMP" ]; then
    echo "ERROR: Missing input file $DUMP - did step 02 (relax) complete?"
    exit 1
fi

# 2. Check python has numpy (no installing - compute nodes are offline)
python3 -c "import numpy; print('python', __import__('sys').version.split()[0], '| numpy', numpy.__version__)"
if [ $? -ne 0 ]; then
    echo "ERROR: python3 on this node cannot import numpy."
    echo "Try 'module avail python' on the login node and add the right"
    echo "'module load' line above, or run this step on your own PC instead:"
    echo "  python 03_region_select_numpy.py $DUMP"
    exit 1
fi

mkdir -p runs/corrosion_stage1

# 3. Run the region extraction
echo "--- Extracting crystal_core / interphase / glass_core patches from $DUMP ---"
python3 03_region_select_numpy.py "$DUMP" 2>&1 | tee runs/corrosion_stage1/log.03_region_select

echo "=== Region Selection Job Finished: $(date) ==="
ls -la slab_*_crystal_core.xyz slab_*_interphase.xyz slab_*_glass_core.xyz 2>&1 || \
    echo "WARNING: expected .xyz files not found - check runs/corrosion_stage1/log.03_region_select"
