#!/bin/bash
#SBATCH --partition=testq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=log.series_model_%j.out
# Geometry-only prediction for every series patch: exposed Ca/Mg/Zn fractions and predicted H2O binding.
module load apps/anaconda3/2024.06
python3 11_adsorption.py analyze > /dev/null
EXTRA=$(ls qe_inputs/c0*_wf.in qe_inputs/c1*_wf.in 2>/dev/null | tr '\n' ' ')
[ -z "$EXTRA" ] && { echo "no series qe_inputs/c*_wf.in found - run prep_series.sh first"; exit 1; }
echo "series cells: $EXTRA"
EXTRA_IN="$EXTRA" python3 12_surface_model.py
