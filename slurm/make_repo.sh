#!/bin/bash
# Assemble the public repository tree.  Run from ~/W3:      bash make_repo.sh
# Writes to $DEST (default ../repo, OUTSIDE W3 so it does not compete with the QE outputs
# for space in the same directory tree).  Set DEST to a scratch filesystem if home is tight:
#     DEST=/scratch/$USER/repo bash make_repo.sh
#
# Errors are NOT suppressed: a failed copy must be visible, not silent.
set -u
DEST="${DEST:-$(cd .. && pwd)/repo}"
NEED_MB=120

echo "=== space"
df -h . | tail -1
AVAIL_MB=$(df -Pm "$(dirname "$DEST")" | tail -1 | awk '{print $4}')
echo "available where the repo will go: ${AVAIL_MB} MB (need about ${NEED_MB} MB)"
if [ "$AVAIL_MB" -lt "$NEED_MB" ]; then
  echo "ERROR: not enough free space. Free some, or set DEST to another filesystem:"
  echo "       DEST=/scratch/\$USER/repo bash make_repo.sh"
  exit 1
fi

rm -rf "$DEST"; mkdir -p "$DEST"/{scripts,slurm,structures/patches,qe_inputs,data/ads,figures} || exit 1
fail=0
grab() {  # grab <destination subdir> <glob...>
  local sub="$1"; shift
  local n=0 b=0
  for f in "$@"; do
    [ -f "$f" ] || continue
    if cp -p "$f" "$DEST/$sub/"; then
      n=$((n+1)); b=$((b + $(stat -c%s "$f")))
    else
      echo "  FAILED to copy $f"; fail=1
    fi
  done
  printf "  %-20s %4d files  %8.1f MB\n" "$sub" "$n" "$(echo "$b" | awk '{print $1/1048576}')"
}

echo "=== copying"
grab scripts [0-9][0-9]_*.py [0-9][0-9]_*.in
grab slurm prep_*.sh run_*.sh submit_*.sh make_repo.sh
grab qe_inputs qe_inputs/*.in
grab structures slab_*.xyz
grab structures/patches patches/*.data
grab data surface_model.csv series_model.csv series.txt analysis_*.csv ads/ads_table.csv
grab data/ads ads/ads_table.csv ads/cref_map.csv ads/sites_*.csv
grab figures figs/*.pdf figs/*.png figs/*.csv

echo "=== verifying (no empty files, nothing over 40 MB)"
EMPTY=$(find "$DEST" -type f -empty | tee /dev/stderr | wc -l)
BIG=$(find "$DEST" -type f -size +40M | tee /dev/stderr | wc -l)
echo "  empty files: $EMPTY        over 40 MB: $BIG"
[ "$EMPTY" -gt 0 ] && { echo "ERROR: empty files above - the copy did not complete. Check space."; fail=1; }

( cd "$DEST" && find . -type f | sort | while read -r f; do
    printf "%10d  %s\n" "$(stat -c%s "$f")" "${f#./}"; done ) > "$DEST/MANIFEST.txt"
( cd "$DEST" && find . -type f ! -name SHA256SUMS -exec sha256sum {} \; | sort -k2 ) > "$DEST/SHA256SUMS"
cat > "$DEST/.gitignore" <<'GI'
*.out
*.err
*.wfc*
*.save/
*.cube
__pycache__/
GI

echo "=== total"
du -sh "$DEST"
if [ "$fail" -eq 0 ]; then
  echo "OK. Repo tree is at: $DEST"
  echo "Add README.md there, then download the folder and push."
else
  echo "FINISHED WITH ERRORS - do not push until they are fixed."
  exit 1
fi
