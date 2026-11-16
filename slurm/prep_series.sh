#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=log.series_%j.out
# Composition series for the merged corrosion paper. NO DFT.
# Cuts glass-matrix patches from each quenched cell, relaxes them with MEAM (box relaxation, which
# 08_relax_patch.in requires for amorphous patches), exports them and writes the QE inputs that
# 12_surface_model.py reads for the geometry-only prediction.
#
#   1. make series.txt: one "<label> <path>" per line, blank lines and # comments ignored
#   2. sbatch prep_series.sh                       (or SERIES_FILE=one.txt sbatch prep_series.sh to test one)
#   3. sbatch submit_series_model.sh
#
# Check in the log that each patch's Mg/Zn/Ca fractions match the composition you expect: that is the
# verification that run_N in W2 is composition N of the Study 2 table.
module purge
module load apps/anaconda3/2024.06
PY="$(command -v python3)"                               # pin it: the lammps module puts its own python first
module load apps/lammps/intel/parallel_studio_xe_2020.4.912/22Jul2025/Optional/1
"$PY" -c "import numpy" 2>/dev/null || { echo "ERROR: $PY has no numpy"; exit 1; }
echo "python: $PY"
SERIES_FILE="${SERIES_FILE:-series.txt}"
# 12 independent, non-overlapping glass patches per composition. Each patch exposes only ~20 atoms within
# 1.5 A of its top, so a single patch is a 20-atom estimate of a surface fraction; pooling 12 gives ~240.
OFFSETS="${OFFSETS:--60,-60 -60,0 -60,60 -20,-60 -20,0 -20,60 20,-60 20,0 20,60 60,-60 60,0 60,60}"
[ -f "$SERIES_FILE" ] || { echo "ERROR: $SERIES_FILE not found"; exit 1; }
LMP="${LMP:-$HOME/lammps_meam/bin/lmp}"                 # same binary submit_patch_relax.sh uses
MPI_BIN="${MPI_BIN:-/opt/apps/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpirun}"
[ -x "$LMP" ] || { echo "ERROR: $LMP not executable. Resubmit with  LMP=/full/path/to/lmp sbatch prep_series.sh"; exit 1; }
for f in library.meam CaZnMg.meam 08_relax_patch.in; do
  [ -f "$f" ] || { echo "ERROR: missing $f in $(pwd)"; exit 1; }
done
NENT=$(grep -cve '^[[:space:]]*#' -e '^[[:space:]]*$' "$SERIES_FILE")
echo "lammps: $LMP"; echo "$SERIES_FILE lists $NENT composition(s)"; mkdir -p patches
NAMES=""
while read -r nm pf _ <&3; do
  case "$nm" in ""|\#*) continue;; esac
  [ -f "$pf" ] || { echo "MISSING parent for $nm: $pf"; continue; }
  for off in $OFFSETS; do
    ox="${off%%,*}"; oy="${off##*,}"
    tag="${nm}_$(echo "$off" | tr ',-' '__')"
    echo; echo "=== $tag   offset $ox $oy   parent $pf"
    "$PY" 07_build_periodic_patch.py build --kind matrix --D 0 --offset "$ox" "$oy" \
            --name "$tag" --parent "$pf" || continue
    "$MPI_BIN" -np $SLURM_NTASKS "$LMP" -in 08_relax_patch.in \
           -var patch_data patches/${tag}.data -var out patches/${tag} -var relax_box 1 \
           -log patches/log.relax_${tag} > /dev/null 2>&1 < /dev/null
    if [ ! -f patches/${tag}_relaxed.data ]; then
      echo "  relax FAILED - last lines of patches/log.relax_${tag}:"; tail -15 patches/log.relax_${tag}; continue
    fi
    grep -E "^FINAL" patches/log.relax_${tag} || echo "  (no FINAL line; check patches/log.relax_${tag})"
    "$PY" 07_build_periodic_patch.py export patches/${tag}_relaxed.data || continue
    NAMES="$NAMES $tag"
  done
done 3< "$SERIES_FILE"     # fd 3: lammps/mpirun read stdin and would otherwise swallow the rest of the list
echo; echo "patches built:$NAMES"
NDONE=$(echo $NAMES | wc -w); NEXP=$((NENT * $(echo $OFFSETS | wc -w)))
[ "$NDONE" -ne "$NEXP" ] && echo "WARNING: $NDONE patches built but $NEXP expected ($NENT compositions x $(echo $OFFSETS | wc -w) offsets)"
if [ -n "$NAMES" ]; then
  "$PY" 05_make_qe_inputs.py $(for n in $NAMES; do echo slab_${n}.xyz; done) 2>&1 | tail -25
  "$PY" 09_make_wf_inputs.py $NAMES 2>&1 | tail -5      # <name>_wf.in is what 12_surface_model.py reads
fi
echo; echo "next:  sbatch submit_series_model.sh"
