#!/bin/bash
# The analysis scripts were written to run in a flat working directory, where ads_table.csv is at
# ads/, the per-atom CSVs are at the top level, and so on. This repository groups those files under
# data/ for readability. Run this once to build the flat layout the scripts expect, using symlinks
# so nothing is duplicated:
#
#     bash setup_workdir.sh && cd work && PYTHONNOUSERSITE=1 python3 ../scripts/22_figures.py make
#
set -eu
R="$(cd "$(dirname "$0")" && pwd)"
W="$R/work"
mkdir -p "$W"/{ads,figs,patches,qe_inputs}
ln -sfn "$R"/data/*.csv "$R"/data/series.txt "$W"/ 2>/dev/null || true
ln -sfn "$R"/data/ads/*.csv "$W"/ads/ 2>/dev/null || true
ln -sfn "$R"/figures/*.csv "$W"/figs/ 2>/dev/null || true
ln -sfn "$R"/structures/*.xyz "$W"/ 2>/dev/null || true
ln -sfn "$R"/structures/patches/*.data "$W"/patches/ 2>/dev/null || true
ln -sfn "$R"/qe_inputs/*.in "$W"/qe_inputs/ 2>/dev/null || true
echo "work/ ready:"
for d in . ads figs patches qe_inputs; do
  printf "  work/%-10s %4d entries\n" "$d" "$(find "$W/$d" -maxdepth 1 -mindepth 1 | wc -l)"
done
echo
echo "Figures need matplotlib built against the same NumPy as the interpreter;"
echo "if a NumPy 2.x user install shadows the system one, prefix with PYTHONNOUSERSITE=1."
