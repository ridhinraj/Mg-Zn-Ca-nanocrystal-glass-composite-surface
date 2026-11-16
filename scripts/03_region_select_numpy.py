"""
Stage 1c (no-OVITO version): extract three DFT patches from the relaxed
equatorial slab using only numpy.

Needs only numpy (no scipy, no OVITO).

Usage:
    python 03_region_select_numpy.py equatorial_slab_relaxed.dump

What it does
  1. Reads the LAMMPS dump (id type x y z ...).
  2. Classifies every atom by its 3D distance r3d from the sphere centre
     (103.5, 103.5, 103.5):  core r<20, interphase 20-30, matrix r>30.
     (3D distance, not lateral: below the cut plane the sphere's
     cross-section shrinks, so lateral distance is only valid AT z=103.5.)
  3. Cuts a square patch (PATCH_SIDE x PATCH_SIDE, DEPTH deep from the
     surface) at three lateral positions along +x from the sphere axis.
  4. Crystallinity cross-check from first-shell order (tight 12-neighbour
     shell + clear gap to the 13th), checked against reference atoms deep in
     the core and deep in the matrix.
  5. Writes one extended-XYZ file per patch (element symbols + lattice),
     ready for 04_xyz_to_qe.py.

Edit the SETTINGS block if a patch comes back far from ~150-250 atoms or
the zone fractions look wrong.
"""

import sys
import numpy as np

# ---------------------------------------------------------------- SETTINGS
COMPOSITION_LABEL = "Mg60Zn35Ca5"
TYPE_TO_ELEMENT = {1: "Mg", 2: "Zn", 3: "Ca"}

CENTER = np.array([103.5, 103.5, 103.5])   # sphere centre
SURFACE_Z = 103.5                          # equatorial cut plane
SURFACE_TOL = 2.0      # A: keep atoms that relaxed UP to this far above the cut plane
                       # (after relaxation the top surface sits ~0.6 A above 103.5)
R_CORE, R_INTERPHASE = 20.0, 30.0          # zone boundaries from bulk profile

PATCH_SIDE = 16.0      # A, square footprint
DEPTH = 16.0           # A below the surface
PAD_XY = 6.0           # A padding added around patch in the exported cell (x,y)
VACUUM_Z = 15.0        # A vacuum above the patch in the exported cell

# lateral offset of each patch centre from the sphere axis, along +x
PATCHES = {
    "crystal_core": 0.0,
    "interphase":   24.0,
    "glass_core":   60.0,
}

# ----------------------------------------------------------------- READ DUMP
def read_dump(fname):
    with open(fname) as f:
        lines = []
        for line in f:
            lines.append(line)
            if line.startswith("ITEM: ATOMS"):
                cols = line.split()[2:]
                break
        n = int(lines[3])
        bx = [list(map(float, lines[5 + i].split()[:2])) for i in range(3)]
        data = np.loadtxt(f)
    ci = {c: i for i, c in enumerate(cols)}
    ids = data[:, ci["id"]].astype(int)
    types = data[:, ci["type"]].astype(int)
    pos = data[:, [ci["x"], ci["y"], ci["z"]]]
    assert len(ids) == n, "atom count mismatch in dump"
    return ids, types, pos, np.array(bx)

# ------------------------------------------------- local-order descriptor
GAP_MIN = 1.15      # crystal: 13th-neighbour distance / 12th-neighbour distance
SPREAD_MAX = 0.06   # crystal: std/mean of the 12 nearest-neighbour distances

def nn_distances(p, L, idx, k=14):
    """Distances to the k nearest neighbours (self excluded) for atoms idx.
    Brute force over all atoms, periodic in x,y only. ~ms per atom, and we
    only query a few thousand atoms, so no neighbour-list library is needed."""
    out = np.empty((len(idx), k - 1))
    for a, i in enumerate(idx):
        d = p - p[i]
        d[:, 0] -= L[0] * np.round(d[:, 0] / L[0])
        d[:, 1] -= L[1] * np.round(d[:, 1] / L[1])
        r2 = np.einsum("ij,ij->i", d, d)
        part = np.argpartition(r2, k)[:k]
        out[a] = np.sqrt(np.sort(r2[part])[1:])
    return out

def is_crystalline(p, L, idx):
    """HCP/FCC-like if the first shell of 12 is tight (low spread) and cleanly
    separated from the second shell (large 13th/12th distance ratio).
    Calibrated on this system: core atoms gap~1.27, spread~0.04;
    matrix atoms gap~1.05, spread~0.07. (Centro-symmetry did NOT separate
    them - ideal HCP has a non-zero CSP baseline ~a^2 - so it is not used.)"""
    D = nn_distances(p, L, idx)
    spread = D[:, :12].std(axis=1) / D[:, :12].mean(axis=1)
    gap = D[:, 12] / D[:, 11]
    return (gap > GAP_MIN) & (spread < SPREAD_MAX)

# ---------------------------------------------------------------------- MAIN
def main(fname):
    ids, types, pos, box = read_dump(fname)
    L = box[:, 1] - box[:, 0]
    print(f"Read {len(ids)} atoms; box x,y side = {L[0]:.2f}, {L[1]:.2f} A")
    print(f"z range of atoms: {pos[:,2].min():.2f} .. {pos[:,2].max():.2f} A")

    # wrap x,y into the box (periodic in x,y; z is a free surface)
    p = pos.copy()
    p[:, 0] = (p[:, 0] - box[0, 0]) % L[0]
    p[:, 1] = (p[:, 1] - box[1, 0]) % L[1]
    p[:, 2] = p[:, 2] - box[2, 0]
    # keep the unwrapped-to-box coordinates for geometry below
    xy_shift = np.array([box[0, 0], box[1, 0], box[2, 0]])
    pw = p + xy_shift            # wrapped positions in original frame

    # periodic-aware r3d from the sphere centre
    d = pw - CENTER
    d[:, 0] -= L[0] * np.round(d[:, 0] / L[0])
    d[:, 1] -= L[1] * np.round(d[:, 1] / L[1])
    r3d = np.linalg.norm(d, axis=1)
    depth = SURFACE_Z - pw[:, 2]

    # ---- sanity-check the crystallinity rule on reference atoms
    rng = np.random.default_rng(0)
    interior = depth > 6.0       # >6 A below the surface -> full neighbour shell
    ref_c = np.where(interior & (r3d < 14.0))[0]
    ref_g = np.where(interior & (r3d > 45.0) & (r3d < 70.0))[0]
    ref_c = rng.choice(ref_c, min(600, len(ref_c)), replace=False)
    ref_g = rng.choice(ref_g, min(600, len(ref_g)), replace=False)
    tpr = 100.0 * is_crystalline(p, L, ref_c).mean()
    fpr = 100.0 * is_crystalline(p, L, ref_g).mean()
    print(f"\nCrystallinity rule check (gap>{GAP_MIN}, spread<{SPREAD_MAX}):")
    print(f"  deep core atoms flagged crystalline   : {tpr:.1f}%  (want high)")
    print(f"  deep matrix atoms flagged crystalline : {fpr:.1f}%  (want low)")
    if tpr < 70 or fpr > 25:
        print("  WARNING: rule separates core/matrix poorly - treat crystal % with caution.")

    # ---- cut patches
    for name, off in PATCHES.items():
        cx = CENTER[0] + off
        cy = CENTER[1]
        dx = d[:, 0] - off
        dy = d[:, 1]
        sel = (np.abs(dx) <= PATCH_SIDE / 2) & (np.abs(dy) <= PATCH_SIDE / 2) \
              & (depth >= -SURFACE_TOL) & (depth <= DEPTH)
        idx = np.where(sel)[0]
        n = len(idx)
        el = np.array([TYPE_TO_ELEMENT[t] for t in types[idx]])
        comp = {e: 100.0 * (el == e).sum() / max(n, 1) for e in ("Mg", "Zn", "Ca")}

        zones = (100.0 * (r3d[idx] < R_CORE).mean(),
                 100.0 * ((r3d[idx] >= R_CORE) & (r3d[idx] <= R_INTERPHASE)).mean(),
                 100.0 * (r3d[idx] > R_INTERPHASE).mean()) if n else (0, 0, 0)

        inner = idx[depth[idx] > 6.0]     # first-shell rule needs a full neighbour shell
        cfrac = 100.0 * is_crystalline(p, L, inner).mean() if len(inner) else float("nan")

        print(f"\n{name}: patch centre offset {off:.0f} A, {n} atoms "
              f"(target ~150-250)")
        print(f"  composition  Mg {comp['Mg']:.1f}%  Zn {comp['Zn']:.1f}%  Ca {comp['Ca']:.1f}%")
        print(f"  by 3D radius core {zones[0]:.0f}%  interphase {zones[1]:.0f}%  matrix {zones[2]:.0f}%")
        print(f"  crystalline among interior atoms: {cfrac:.1f}%  "
              f"({len(inner)} atoms tested)")

        # ---- export extended xyz with a periodic cell around the patch
        # positions relative to the patch (handles periodic wrap near edges)
        rel = np.column_stack([dx[idx], dy[idx], pw[idx, 2]])
        rel[:, 0] -= rel[:, 0].min(); rel[:, 1] -= rel[:, 1].min()
        rel[:, 2] -= rel[:, 2].min()
        cell = np.array([rel[:, 0].max() + 2 * PAD_XY,
                         rel[:, 1].max() + 2 * PAD_XY,
                         rel[:, 2].max() + VACUUM_Z])
        rel[:, 0] += PAD_XY; rel[:, 1] += PAD_XY
        fn = f"slab_{COMPOSITION_LABEL}_{name}.xyz"
        with open(fn, "w") as f:
            f.write(f"{n}\n")
            f.write(f'Lattice="{cell[0]:.4f} 0 0 0 {cell[1]:.4f} 0 0 0 {cell[2]:.4f}" '
                    f'Properties=species:S:1:pos:R:3 pbc="T T T"\n')
            for e, r in zip(el, rel):
                f.write(f"{e} {r[0]:.5f} {r[1]:.5f} {r[2]:.5f}\n")
        print(f"  wrote {fn}  (cell {cell[0]:.1f} x {cell[1]:.1f} x {cell[2]:.1f} A)")

    print("\nDone. Check atom counts, zone fractions and crystal % before "
          "moving on to 04_xyz_to_qe.py.")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "equatorial_slab_relaxed.dump")
