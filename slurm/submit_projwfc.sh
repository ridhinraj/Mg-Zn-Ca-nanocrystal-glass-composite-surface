#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=40G
#SBATCH --time=01:00:00
#SBATCH --output=log.projwfc_%j.out
#SBATCH --error=log.projwfc_%j.err

# Projected DOS + Lowdin charges from a finished pw.x SCF (needs tmp_<NAME>/<NAME>.save).
# Light job for the Gamma-only crystal/interface cells (minutes). It is NOT light for the big cells:
#  - 3x3x1 k-point runs (5 k-points) and the matrix patches (~950 atomic wavefunctions, ~830 bands)
#    can take well over 30 minutes and several GB. testq allows at most 1 hour, so for those use e.g.
#      sbatch --partition=mediumq --time=04:00:00 --mem=60G --export=ALL,NAME=matrix_xp submit_projwfc.sh
# Patch is chosen at submit time:
#   sbatch --export=ALL,NAME=Mg60Zn35Ca5_interphase submit_projwfc.sh
NAME="${NAME:-Mg60Zn35Ca5_crystal_core}"
SAVE="tmp_${NAME}/${NAME}.save"

echo "=== projwfc started: $(date) | patch: $NAME ==="
module purge
module load apps/quantum-espresso/7.5/cpu
export OMP_NUM_THREADS=1

if [ ! -d "$SAVE" ]; then
    echo "ERROR: $SAVE not found - did pw.x finish with 'Writing all to output data dir'?"; exit 1
fi
mkdir -p pdos/${NAME}

# Energy window follows this run's Fermi energy (E_F differs between cells: about -2 eV for the
# island patches but about +0.3 eV for the periodic ones). Energies in eV, absolute. degauss in Ry.
EF=$(grep "the Fermi energy is" qe_inputs/${NAME}.out | tail -1 | awk '{print $5}')
EMIN=$(awk -v e=$EF 'BEGIN{printf "%.2f", e-12.0}')
EMAX=$(awk -v e=$EF 'BEGIN{printf "%.2f", e+4.0}')
echo "E_F = $EF eV -> PDOS window $EMIN .. $EMAX eV"
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

mpirun -np $SLURM_NTASKS projwfc.x -nk 1 -in qe_inputs/${NAME}.projwfc.in > qe_inputs/${NAME}.projwfc.out 2>&1

echo "=== projwfc finished: $(date) ==="
grep -E "JOB DONE" qe_inputs/${NAME}.projwfc.out || echo "WARNING: no JOB DONE - see qe_inputs/${NAME}.projwfc.out"
ls pdos/${NAME} | head -5; echo "files in pdos/${NAME}: $(ls pdos/${NAME} | wc -l)"
