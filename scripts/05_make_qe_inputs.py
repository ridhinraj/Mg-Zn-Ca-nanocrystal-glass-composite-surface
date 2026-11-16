"""
Stage 2a: build Quantum ESPRESSO pw.x inputs from the patch .xyz files
written by 03_region_select_numpy.py.  Plain Python, no ASE / numpy needed.

Usage (from W3/):
    python3 05_make_qe_inputs.py                  # all slab_*.xyz in this folder
    python3 05_make_qe_inputs.py slab_Mg60Zn35Ca5_crystal_core.xyz

Writes qe_inputs/<name>.in and prints, for each patch, the composition, the
valence-electron count (read from the UPF headers) and the default band count.

Design choices (edit the SETTINGS block to change them):
  * calculation = scf on the MEAM-relaxed geometry (no DFT relaxation).
    COMPUTE_FORCES (off by default) prints DFT forces, which show how far the MEAM
    geometry is from DFT equilibrium - worth reporting in the thesis later.
  * Gamma-point only ('K_POINTS gamma'): defensible for a disordered finite
    patch, and QE's gamma trick roughly halves cost and memory.
  * Marzari-Vanderbilt smearing 0.02 Ry (the SSSP-efficiency convention).
  * mixing_mode local-TF: recommended for inhomogeneous / slab systems.
  * Cutoffs are a STARTING POINT: 40 / 200 Ry, the values GBRV potentials were
    designed for (Mg PAW header asks for >= 87 Ry density cutoff). Confirm with a
    cutoff test (e.g. raise ecutrho to 320) before production runs.
"""

import glob
import os
import re
import sys

# ---------------------------------------------------------------- SETTINGS
PSEUDO_DIR = os.environ.get("PSEUDO_DIR",
                            "/gpfs-home/p230674mt/W3/SSSP_1.3.0_PBE_efficiency")
PSEUDOS = {
    "Mg": "Mg.pbe-n-kjpaw_psl.0.3.0.UPF",
    "Zn": "Zn_pbe_v1.uspp.F.UPF",
    "Ca": "Ca_pbe_v1.uspp.F.UPF",
}
MASSES = {"Mg": 24.305, "Zn": 65.38, "Ca": 40.078}

ECUTWFC = 40.0          # Ry   (starting point - confirm with a cutoff test)
ECUTRHO = 200.0         # Ry   (GBRV design value, dual 5; 320 gave a 288x288x320 grid and ~63 GB)
DEGAUSS = 0.02          # Ry, Marzari-Vanderbilt smearing
CONV_THR = 1.0e-6       # Ry total (~6e-9 per atom; plenty for PDOS / charges)
MIXING_BETA = 0.15      # 0.3 gave charge sloshing (accuracy jumped 0.6 -> 19 Ry, 47 iterations)
COMPUTE_FORCES = False  # QE writes the saved wavefunctions/density only AFTER the forces step,
                        # so a crash there loses the whole SCF. Off for now; the forces are a
                        # nice-to-have (MEAM-vs-DFT check) and can be run later from the saved data.
DIPOLE_CORRECTION = False   # slab dipole correction (tefield/dipfield); untested here
VACUUM_Z = 15.0         # A, only used to place the dipole sawtooth
# k-points: default Gamma only. For the k-point test on periodic cells:
#   KPTS="2 2 1" TAG=k221 python3 05_make_qe_inputs.py slab_crystal_ideal.xyz
# writes qe_inputs/crystal_ideal_k221.in (so runs do not overwrite each other).
# Extra bands above the occupied ones (PDOS reach above E_F). ~50 bands per eV near E_F for these
# 16x17 A cells, so 110 extra reaches about E_F + 2 eV. Override: NBND_EXTRA=160 for ~ +3 eV, 0 = QE default.
NBND_EXTRA = int(os.environ.get("NBND_EXTRA", "110"))
KPTS = os.environ.get("KPTS", "").split()
TAG = os.environ.get("TAG", "")


# ------------------------------------------------------------------ helpers
def read_zval(upf_path):
    """Valence charge from a UPF header (handles UPF v2 xml and old v1)."""
    try:
        text = open(upf_path, errors="ignore").read(20000)
    except OSError:
        return None
    m = re.search(r'z_valence\s*=\s*"\s*([0-9.eEdD+-]+)\s*"', text)
    if m:
        return float(m.group(1).replace("D", "E").replace("d", "e"))
    m = re.search(r"([0-9.eEdD+-]+)\s+Z valence", text)
    if m:
        return float(m.group(1).replace("D", "E").replace("d", "e"))
    return None


def read_xyz(fname):
    with open(fname) as f:
        n = int(f.readline())
        header = f.readline()
        m = re.search(r'Lattice="([^"]+)"', header)
        if not m:
            raise ValueError(f"{fname}: no Lattice= in the comment line")
        lat = [float(x) for x in m.group(1).split()]
        cell = [lat[0:3], lat[3:6], lat[6:9]]
        atoms = []
        for _ in range(n):
            s, x, y, z = f.readline().split()[:4]
            atoms.append((s, float(x), float(y), float(z)))
    return cell, atoms


def build_input(name, cell, atoms, nbnd=None):
    species = sorted({a[0] for a in atoms}, key=lambda e: list(PSEUDOS).index(e))
    nat = len(atoms)

    system_extra = ""
    if DIPOLE_CORRECTION:
        cz = cell[2][2]
        ztop = max(a[3] for a in atoms)
        emaxpos = min((ztop + 0.5 * VACUUM_Z) / cz, 0.98)
        system_extra = (f"  tefield = .true.\n  dipfield = .true.\n  edir = 3\n"
                        f"  emaxpos = {emaxpos:.3f}\n  eopreg = 0.08\n  eamp = 0.0\n")

    out = []
    out.append("&control\n  calculation = 'scf'\n"
               f"  prefix = '{name}'\n  outdir = './tmp_{name}'\n"
               f"  pseudo_dir = '{PSEUDO_DIR}'\n"
               "  tprnfor = %s\n  disk_io = 'low'\n/\n" % (".true." if COMPUTE_FORCES else ".false."))
    out.append("&system\n  ibrav = 0\n"
               f"  nat = {nat}\n  ntyp = {len(species)}\n"
               f"  ecutwfc = {ECUTWFC}\n  ecutrho = {ECUTRHO}\n"
               + (f"  nbnd = {nbnd}\n" if nbnd else "") +
               "  occupations = 'smearing'\n  smearing = 'mv'\n"
               f"  degauss = {DEGAUSS}\n{system_extra}/\n")
    out.append("&electrons\n  conv_thr = %.1e\n  mixing_beta = %s\n"
               "  mixing_mode = 'local-TF'\n  electron_maxstep = 200\n/\n"
               % (CONV_THR, MIXING_BETA))
    out.append("ATOMIC_SPECIES\n" + "".join(
        f"  {e}  {MASSES[e]}  {PSEUDOS[e]}\n" for e in species))
    out.append("CELL_PARAMETERS angstrom\n" + "".join(
        "  %.6f  %.6f  %.6f\n" % tuple(v) for v in cell))
    out.append("ATOMIC_POSITIONS angstrom\n" + "".join(
        "  %s  %.6f  %.6f  %.6f\n" % a for a in atoms))
    if KPTS:
        out.append("K_POINTS automatic\n  %s %s %s 0 0 0\n" % tuple(KPTS))
    else:
        out.append("K_POINTS gamma\n")
    return "\n".join(out)


# ---------------------------------------------------------------------- main
def main(files):
    os.makedirs("qe_inputs", exist_ok=True)
    zval = {e: read_zval(os.path.join(PSEUDO_DIR, p)) for e, p in PSEUDOS.items()}
    print("Valence electrons per atom read from the UPF headers:", zval)
    if any(v is None for v in zval.values()):
        print("  WARNING: could not read z_valence for some elements "
              "(is PSEUDO_DIR correct?). Electron counts below will be incomplete.")

    for fn in files:
        name = os.path.basename(fn).replace(".xyz", "").replace("slab_", "") + (f"_{TAG}" if TAG else "")
        cell, atoms = read_xyz(fn)
        counts = {}
        for a in atoms:
            counts[a[0]] = counts.get(a[0], 0) + 1
        missing = [e for e in counts if e not in PSEUDOS]
        if missing:
            raise SystemExit(f"{fn}: no pseudopotential defined for {missing}")
        nbnd = None
        if NBND_EXTRA > 0 and all(zval[e] is not None for e in counts):
            nel_ = sum(zval[e] * c for e, c in counts.items())
            nocc_ = int(round(nel_ / 2))
            qe_default = max(round(1.2 * nocc_), nocc_ + 4)       # what pw.x picks if nbnd is not given
            nbnd = max(qe_default, nocc_ + NBND_EXTRA)           # never FEWER than the default
        with open(f"qe_inputs/{name}.in", "w") as f:
            f.write(build_input(name, cell, atoms, nbnd))

        msg = f"{name}: {len(atoms)} atoms {counts}, cell {cell[0][0]:.1f} x {cell[1][1]:.1f} x {cell[2][2]:.1f} A"
        if all(zval[e] is not None for e in counts):
            nel = sum(zval[e] * c for e, c in counts.items())
            nbnd = max(round(1.2 * nel / 2), round(nel / 2) + 4)   # QE default, smearing
            nocc = int(round(nel / 2)); dflt = max(round(1.2 * nocc), nocc + 4)
            nb = max(dflt, nocc + NBND_EXTRA) if NBND_EXTRA > 0 else dflt
            msg += f"\n   -> {nel:.0f} valence electrons, nbnd = {nb} (occupied {nocc} + {nb - nocc} extra; QE default would be {dflt})"
        print(msg)
        print(f"   wrote qe_inputs/{name}.in")


if __name__ == "__main__":
    main(sys.argv[1:] or sorted(glob.glob("slab_*.xyz")))
