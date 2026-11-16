"""
Stage 2 (Option B): laterally PERIODIC surface patches for DFT.

Why: finite "island" patches have under-coordinated edge atoms whose electronic
descriptors are shifted by 0.4-0.5 eV (more than the surface effect itself) and
only ~30% of the atoms are usable. A periodic x,y cell has no edges.

Cell: the nanocrystal measured in AN_composite_relaxed_300K.data is HCP with
a = 3.245 A, c = 5.300 A, c-axis along z, a-axis along x. A periodic cell is
seamless for the crystal only if its sides are whole multiples of the lattice
period, so  Lx = 5a = 16.225 A,  Ly = 3*sqrt(3)*a = 16.862 A  (30 atoms/plane).

Three kinds of patch
  crystal   ideal Mg(0001) slab, 6 planes (no parent needed).
  interphase  cut from the parent with the surface plane D angstrom from the
            sphere centre, cell centred on the polar axis. On the polar axis the
            crystal is pure Mg up to D ~ 24 A and the matrix starts at D ~ 28 A, so
            D = 24..28 is a "corrosion front approaching the nanocrystal" series.
            --pole south uses the opposite pole (180-degree rotation about x,
            keeps the lattice orientation): independent glass atoms, equivalent
            crystal termination -> a replicate.
  matrix    cut at the sphere's equatorial height, cell far from the axis
            (--offset 60 0): amorphous matrix only. Offsets +-60 A in x and y give
            4 independent (non-overlapping) replicates.

Usage
  python3 07_build_periodic_patch.py series  --parent AN_composite_relaxed_300K.data
  python3 07_build_periodic_patch.py build --kind interphase --D 26 --pole north \\
          --name iface_N_D26 --parent AN_composite_relaxed_300K.data
  python3 07_build_periodic_patch.py export patches/iface_N_D26_relaxed.data
Outputs go to patches/: <name>.data (LAMMPS, to relax with 08_relax_patch.in) and
slab_<name>.xyz (QE-ready after relaxation: use `export`).
"""

import argparse
import os
import sys

import numpy as np

# ---------------------------------------------------------------- SETTINGS
A_LAT, C_LAT = 3.245, 5.300          # measured in the parent nanocrystal
NX, NY = 5, 3
LX, LY = NX * A_LAT, NY * np.sqrt(3.0) * A_LAT
THICKNESS = 14.0                     # A of material below the surface plane
VACUUM_Z = 15.0                      # A of vacuum above the top surface (xyz export)
LMP_ZPAD_BOTTOM, LMP_ZPAD_TOP = 3.0, 20.0   # box padding in the LAMMPS data file
MIN_DIST = 2.2                       # A, atoms closer than this (periodic) are overlaps
CENTER = (103.5, 103.5, 103.5)       # nanocrystal centre in the parent
MASSES = {1: 24.32, 2: 65.39, 3: 40.078}
EL = {1: "Mg", 2: "Zn", 3: "Ca"}
SYM_TO_TYPE = {v: k for k, v in EL.items()}


# ------------------------------------------------------------------ I/O
def read_parent(fn):
    with open(fn) as f:
        lines = f.readlines()
    n = int(lines[2].split()[0])
    lo_hi = []
    for l in lines[:20]:
        if "xlo" in l or "ylo" in l or "zlo" in l:
            lo_hi.append([float(v) for v in l.split()[:2]])
    s = next(i for i, l in enumerate(lines) if l.strip().startswith("Atoms")) + 2
    arr = np.array([l.split()[:5] for l in lines[s:s + n]], dtype=float)
    return arr[:, 1].astype(int), arr[:, 2:5], np.array(lo_hi)


def write_lammps(path, types, pos, lx, ly, zlo, zhi):
    with open(path, "w") as f:
        f.write(f"periodic patch ({len(types)} atoms)\n\n{len(types)} atoms\n3 atom types\n\n")
        f.write(f"0.0 {lx:.6f} xlo xhi\n0.0 {ly:.6f} ylo yhi\n{zlo:.6f} {zhi:.6f} zlo zhi\n\n")
        f.write("Masses\n\n" + "".join(f"{t} {m}\n" for t, m in MASSES.items()) + "\nAtoms # atomic\n\n")
        for i, (t, p) in enumerate(zip(types, pos), 1):
            f.write(f"{i} {t} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")


def write_xyz(path, types, pos, lx, ly):
    pos = pos.copy()
    pos[:, 0] %= lx
    pos[:, 1] %= ly
    pos[:, 2] -= pos[:, 2].min()
    lz = pos[:, 2].max() + VACUUM_Z
    with open(path, "w") as f:
        f.write(f"{len(types)}\n")
        f.write(f'Lattice="{lx:.4f} 0 0 0 {ly:.4f} 0 0 0 {lz:.4f}" '
                f'Properties=species:S:1:pos:R:3 pbc="T T T"\n')
        for t, p in zip(types, pos):
            f.write(f"{EL[t]} {p[0]:.5f} {p[1]:.5f} {p[2]:.5f}\n")


# ---------------------------------------------------------------- geometry
def min_image_xy(d, lx, ly):
    d = d.copy()
    d[..., 0] -= lx * np.round(d[..., 0] / lx)
    d[..., 1] -= ly * np.round(d[..., 1] / ly)
    return d


def remove_overlaps(types, pos, lx, ly):
    """Iteratively delete atoms involved in contacts shorter than MIN_DIST (periodic x,y)."""
    removed = 0
    while True:
        d = min_image_xy(pos[:, None, :] - pos[None, :, :], lx, ly)
        r = np.linalg.norm(d, axis=2)
        np.fill_diagonal(r, 1e9)
        close = r < MIN_DIST
        if not close.any():
            return types, pos, removed
        worst = np.argmax(close.sum(axis=1))
        keep = np.ones(len(types), bool)
        keep[worst] = False
        types, pos = types[keep], pos[keep]
        removed += 1


def build_crystal():
    """Ideal HCP Mg(0001) slab, 6 planes, in the commensurate NX x NY orthorhombic cell."""
    a, c = A_LAT, C_LAT
    A = [(0.0, 0.0), (a / 2, np.sqrt(3) * a / 2)]
    cells = [(i * a, j * np.sqrt(3) * a) for i in range(NX) for j in range(NY)]
    pos = []
    for k in range(6):
        shift = (0.0, 0.0) if k % 2 == 0 else (0.0, a / np.sqrt(3))
        for cx, cy in cells:
            for ax, ay in A:
                pos.append(((cx + ax + shift[0]) % LX, (cy + ay + shift[1]) % LY, k * c / 2))
    pos = np.array(pos)
    return np.ones(len(pos), int), pos


def build_from_parent(parent, kind, D, pole, offset, thickness, center):
    types, pos, box = parent
    L = box[:, 1] - box[:, 0]
    d = pos - np.array(center)
    d -= L * np.round(d / L)                      # minimum image about the centre
    if pole == "south":                           # 180 deg rotation about x: (x,y,z)->(x,-y,-z)
        d[:, 1] *= -1
        d[:, 2] *= -1
    ox, oy = offset
    dx, dy = d[:, 0] - ox, d[:, 1] - oy
    sel = ((dx >= -LX / 2) & (dx < LX / 2) & (dy >= -LY / 2) & (dy < LY / 2)
           & (d[:, 2] <= D) & (d[:, 2] > D - thickness))
    t = types[sel]
    p = np.column_stack([dx[sel] + LX / 2, dy[sel] + LY / 2, d[sel, 2] - (D - thickness)])
    return t, p


def summarize(label, t, p, removed=None):
    f = [100.0 * (t == k).mean() for k in (1, 2, 3)]
    top = p[:, 2] > p[:, 2].max() - 3.0
    ft = [100.0 * (t[top] == k).mean() for k in (1, 2, 3)]
    nel = 2 * (t == 1).sum() + 20 * (t == 2).sum() + 10 * (t == 3).sum()
    s = (f"{label:14s} {len(t):4d} atoms  Mg/Zn/Ca = {f[0]:4.1f}/{f[1]:4.1f}/{f[2]:4.1f}  "
         f"top 3 A: {ft[0]:3.0f}/{ft[1]:3.0f}/{ft[2]:3.0f}  ~{nel} valence e-")
    if removed:
        s += f"  (removed {removed} overlapping)"
    print(s)


def make(name, kind, D, pole, offset, parent, thickness, center, outdir="patches"):
    os.makedirs(outdir, exist_ok=True)
    if kind == "crystal":
        t, p = build_crystal()
        removed = 0
    else:
        t, p = build_from_parent(parent, kind, D, pole, offset, thickness, center)
        t, p, removed = remove_overlaps(t, p, LX, LY)
    zmax = p[:, 2].max()
    write_lammps(f"{outdir}/{name}.data", t, p, LX, LY, -LMP_ZPAD_BOTTOM, zmax + LMP_ZPAD_TOP)
    write_xyz(f"{outdir}/slab_{name}.xyz", t, p, LX, LY)      # UNrelaxed, for inspection
    summarize(name, t, p, removed)
    return t, p


# ---------------------------------------------------------------- commands
def cmd_export(path):
    """Relaxed LAMMPS data file -> slab_<name>.xyz (wrapped, ready for 05_make_qe_inputs.py)."""
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[2].split()[0])
    lx = ly = None
    for l in lines[:20]:
        if "xlo" in l:
            lx = float(l.split()[1]) - float(l.split()[0])
        if "ylo" in l:
            ly = float(l.split()[1]) - float(l.split()[0])
    s = next(i for i, l in enumerate(lines) if l.strip().startswith("Atoms")) + 2
    arr = np.array([l.split()[:5] for l in lines[s:s + n]], dtype=float)
    o = np.argsort(arr[:, 0])
    arr = arr[o]
    name = os.path.basename(path).replace("_relaxed.data", "").replace(".data", "")
    out = f"slab_{name}.xyz"
    t, p = arr[:, 1].astype(int), arr[:, 2:5]
    write_xyz(out, t, p, lx, ly)
    # quality checks on the relaxed patch (periodic contacts, density, composition)
    d = p[:, None, :] - p[None, :, :]
    d[..., 0] -= lx * np.round(d[..., 0] / lx)
    d[..., 1] -= ly * np.round(d[..., 1] / ly)
    r = np.linalg.norm(d, axis=2)
    np.fill_diagonal(r, 1e9)
    thick = p[:, 2].max() - p[:, 2].min()
    comp = [100.0 * (t == k).mean() for k in (1, 2, 3)]
    print(f"wrote {out}  ({n} atoms, cell {lx:.3f} x {ly:.3f} A, thickness {thick:.1f} A)")
    print(f"  composition Mg/Zn/Ca = {comp[0]:.1f}/{comp[1]:.1f}/{comp[2]:.1f}   "
          f"density {n / (lx * ly * thick):.4f} atoms/A^3   shortest contact {r.min():.2f} A")
    if r.min() < 2.3:
        print("  WARNING: contact < 2.3 A - relaxation did not resolve an overlap; check the log")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--kind", choices=["crystal", "interphase", "matrix"], required=True)
    b.add_argument("--name", required=True)
    b.add_argument("--parent", default="AN_composite_relaxed_300K.data")
    b.add_argument("--D", type=float, default=26.0)
    b.add_argument("--pole", choices=["north", "south"], default="north")
    b.add_argument("--offset", type=float, nargs=2, default=[0.0, 0.0])
    b.add_argument("--thickness", type=float, default=THICKNESS)
    s = sub.add_parser("series")
    s.add_argument("--parent", default="AN_composite_relaxed_300K.data")
    e = sub.add_parser("export")
    e.add_argument("data")
    args = ap.parse_args()

    if args.cmd == "export":
        return cmd_export(args.data)
    print(f"cell {LX:.3f} x {LY:.3f} A (a={A_LAT}, c={C_LAT}); thickness {THICKNESS} A\n")
    if args.cmd == "build":
        parent = read_parent(args.parent) if args.kind != "crystal" else None
        make(args.name, args.kind, args.D, args.pole, args.offset, parent, args.thickness, CENTER)
    else:
        parent = read_parent(args.parent)
        make("crystal_ideal", "crystal", 0, "north", (0, 0), None, THICKNESS, CENTER)
        for pole, tag in (("north", "N"), ("south", "S")):
            for D in (24, 25, 26, 27, 28):
                make(f"iface_{tag}_D{D}", "interphase", D, pole, (0, 0), parent, THICKNESS, CENTER)
        for tag, off in (("xp", (60, 0)), ("xm", (-60, 0)), ("yp", (0, 60)), ("ym", (0, -60))):
            make(f"matrix_{tag}", "matrix", 0.0, "north", off, parent, THICKNESS, CENTER)
        print("\nRecommended pilot: crystal_ideal, iface_N_D26, iface_S_D26, matrix_xp, matrix_xm")


if __name__ == "__main__":
    main()
