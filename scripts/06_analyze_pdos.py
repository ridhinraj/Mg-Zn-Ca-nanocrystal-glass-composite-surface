"""
Stage 2b: per-atom electronic descriptors from a finished pw.x + projwfc.x run.

Usage (from W3/, plain Python + numpy; fine on the login node or your PC):
    python3 06_analyze_pdos.py Mg60Zn35Ca5_crystal_core [geometry.xyz]

Needs:
    slab_<NAME>.xyz                 patch geometry (atom order = QE atom order)
    qe_inputs/<NAME>.out            pw.x output (for the Fermi energy)
    qe_inputs/<NAME>.projwfc.out    projwfc.x output (Lowdin charges)
    pdos/<NAME>/<NAME>.pdos_atm#N(El)_wfc#M(l)    projected DOS files

Per atom it computes
    q            Lowdin charge (electrons) projected on the atomic s,p,(d) states
    dq           q minus the mean q of the same element in the deep interior
                 (> INTERIOR_DEPTH below the surface) of the SAME patch.
                 Lowdin totals are NOT directly comparable to the nominal valence
                 (projections miss some charge), so charge is reported relative
                 to interior atoms, not to the valence number.
    N_s, N_p     local DOS at the Fermi level (states/eV), averaged over +-EF_HALF_WIN
    e_p_occ      p-band centre of the OCCUPIED part: first moment of p-PDOS from the
                 band bottom to E_F, relative to E_F (eV)
    e_p_win      p-band centre over [E_F+WIN_LO, E_F+WIN_HI] (eV, relative to E_F)
and prints mean +- std grouped by element and by depth below the top surface.
Everything is also written to analysis_<NAME>.csv.

CAUTION: the upper edge of any window is limited by how many bands the SCF
included. Check the "highest PDOS energy" line printed below before trusting
e_p_win; e_p_occ and N(E_F) do not have this limitation.
"""

import glob
import os
import re
import sys

import numpy as np

# ---------------------------------------------------------------- SETTINGS
EF_HALF_WIN = 0.25          # eV, window for DOS at E_F
WIN_LO, WIN_HI = -3.0, 1.5  # eV relative to E_F, for e_p_win (needs bands up to ~E_F+2 eV: NBND_EXTRA=110 in 05)
SURFACE_DEPTH = float(os.environ.get("SURFACE_DEPTH", "3.0"))  # A below the top atom: "surface" = outermost ~2 atomic planes.
                            # In amorphous/Zn-capped cells the first Mg sits 2-3 A below the topmost (Zn) atom, so 2.0 missed it.
SUBSURFACE_DEPTH = 6.0      # A: sub-surface ends here; beyond = "interior"
INTERIOR_DEPTH = 6.0        # A: reference atoms for dq
BOTTOM_MARGIN = 2.0         # A above the lowest atom: the patch's BOTTOM face is a free surface made by the
                            # cut, not a real interface -> excluded from the statistics (kept in the CSV, bottom=1).
NN_CUT = 3.4               # A: a Zn/Ca atom within this distance counts as a "solute neighbour" (n_ZnCa)
ZVAL = {"Mg": 2.0, "Zn": 20.0, "Ca": 10.0}   # nominal valence of the SSSP pseudopotentials (dZ = q - ZVAL)
EDGE_MARGIN = 3.0           # A: atoms this close to the patch's lateral edge are excluded from
                            # the summary statistics (finite island -> under-coordinated). They stay
                            # in the CSV with edge=1.


# ------------------------------------------------------------------ parsers
def parse_ef(pw_out):
    for line in open(pw_out, errors="ignore"):
        m = re.search(r"the Fermi energy is\s+([-0-9.]+)\s*ev", line, re.I)
        if m:
            ef = float(m.group(1))
    return ef


def parse_top_band(pw_out):
    """Highest Kohn-Sham eigenvalue (eV) printed after 'End of self-consistent calculation'."""
    vals, on = [], False
    for line in open(pw_out, errors="ignore"):
        if "bands (ev)" in line:
            on = True
            continue
        if on and "the Fermi energy" in line:
            break
        if on:
            vals += [float(x) for x in re.findall(r"-?\d+\.\d+", line)]
    return max(vals) if vals else None


def parse_lowdin(path):
    q = {}
    pat = re.compile(r"Atom #\s*(\d+): total charge =\s*([-0-9.]+),\s*([spdf])\s*=\s*([-0-9.]+)")
    started = False
    for line in open(path, errors="ignore"):
        if "Lowdin Charges" in line:
            started = True
        if not started:
            continue
        m = pat.search(line)
        if m:
            i = int(m.group(1))
            d = q.setdefault(i, {"tot": float(m.group(2))})
            d[m.group(3)] = float(m.group(4))
    return q


def read_xyz(path):
    with open(path) as f:
        n = int(f.readline())
        header = f.readline()
        m = re.search(r'Lattice="([^"]+)"', header)
        lat = [float(x) for x in m.group(1).split()] if m else None
        el, pos = [], []
        for _ in range(n):
            s, x, y, z = f.readline().split()[:4]
            el.append(s)
            pos.append([float(x), float(y), float(z)])
    return el, np.array(pos), lat


def load_pdos(folder, name):
    pat = re.compile(r"pdos_atm#(\d+)\((\w+)\)_wfc#(\d+)\((\w)\)")
    data, first = {}, None
    for fn in glob.glob(os.path.join(folder, f"{name}.pdos_atm#*")):
        m = pat.search(os.path.basename(fn))
        if not m:
            continue
        a, l = int(m.group(1)), m.group(4)
        if first is None:
            first = fn
            print("First PDOS file:", os.path.basename(fn))
            print("  header:", open(fn).readline().strip())
        arr = np.loadtxt(fn, comments="#")
        # Zn and Ca have SEVERAL shells of the same l (semicore 3s/3p plus valence 4s/4p, and 3d): one file each.
        # Sum all files of the same l (column 1 = m-summed LDOS). Semicore shells lie far below the PDOS window
        # (E_F-12 eV), so they add nothing there; before this fix the last file read silently overwrote the others.
        if l in data.setdefault(a, {}):
            E0, r0 = data[a][l]
            data[a][l] = (E0, r0 + arr[:, 1])
        else:
            data[a][l] = (arr[:, 0], arr[:, 1].copy())     # E (eV), ldos(E)
    return data


# ---------------------------------------------------------------- numerics
def band_centre(E, rho, lo, hi):
    m = (E >= lo) & (E <= hi)
    w = rho[m].sum()
    return (E[m] * rho[m]).sum() / w if w > 1e-8 else np.nan


def main(name, xyz=None):
    ef = parse_ef(f"qe_inputs/{name}.out")
    import os, warnings
    warnings.filterwarnings("ignore", message="Mean of empty slice")
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    pf = f"qe_inputs/{name}.projwfc.out"
    if not os.path.exists(pf):
        sys.exit(f"ERROR: {pf} does not exist - projwfc has not run for NAME={name}. "
                 f"Check: squeue -u $USER ; sacct -u $USER --format=JobID,JobName%22,State,Elapsed,MaxRSS | grep -i projwfc")
    q = parse_lowdin(pf)
    el, pos, lat = read_xyz(xyz or f"slab_{name}.xyz")
    # periodic cell (atoms span nearly the whole x,y cell) -> no lateral edge to exclude
    periodic = lat is not None and (lat[0] - np.ptp(pos[:, 0]) < 3.0 or lat[4] - np.ptp(pos[:, 1]) < 3.0) \
        and (lat[0] < 40.0)
    edge_margin = 0.0 if periodic else EDGE_MARGIN
    print("periodic x,y cell: lateral-edge exclusion OFF" if periodic else "finite island patch: lateral-edge exclusion ON")
    pdos = load_pdos(f"pdos/{name}", name)
    n = len(el)
    if not pdos or not q:
        sys.exit(f"ERROR: no PDOS files in pdos/{name}/ or no Lowdin block in qe_inputs/{name}.projwfc.out - "
                 f"has submit_projwfc.sh finished for NAME={name}? (check squeue and the projwfc output)")
    print(f"\n{name}: {n} atoms, E_F = {ef} eV, {len(q)} Lowdin entries, {len(pdos)} atoms with PDOS")
    if len(q) != n or len(pdos) != n:
        print("  WARNING: atom counts disagree between xyz / Lowdin / PDOS - check the inputs")

    depth = pos[:, 2].max() - pos[:, 2]
    x0, x1, y0, y1 = pos[:, 0].min(), pos[:, 0].max(), pos[:, 1].min(), pos[:, 1].max()
    edge_dist = np.minimum.reduce([pos[:, 0] - x0, x1 - pos[:, 0], pos[:, 1] - y0, y1 - pos[:, 1]])
    # solute-neighbour count: Zn/Ca atoms within NN_CUT (periodic in x,y for periodic cells)
    is_sol = np.array([e != "Mg" for e in el])
    nn = np.zeros(n, int)
    for i in range(n):
        dv = pos[is_sol] - pos[i]
        if periodic:
            dv[:, 0] -= lat[0] * np.round(dv[:, 0] / lat[0])
            dv[:, 1] -= lat[4] * np.round(dv[:, 1] / lat[4])
        dist = np.sqrt((dv ** 2).sum(1))
        nn[i] = int(((dist < NN_CUT) & (dist > 0.1)).sum())
    rows = []
    for i in range(1, n + 1):
        d = pdos.get(i, {})
        r = {"atom": i, "el": el[i - 1], "x": pos[i - 1, 0], "y": pos[i - 1, 1],
             "z": pos[i - 1, 2], "depth": depth[i - 1], "edge": int(edge_dist[i - 1] < edge_margin), "bottom": int(pos[i - 1, 2] - pos[:, 2].min() < BOTTOM_MARGIN),
             "n_ZnCa": int(nn[i - 1]), "q": q.get(i, {}).get("tot", np.nan),
             "q_s": q.get(i, {}).get("s", np.nan), "q_p": q.get(i, {}).get("p", np.nan)}
        for l in ("s", "p"):
            if l in d:
                E, rho = d[l]
                m = (E >= ef - EF_HALF_WIN) & (E <= ef + EF_HALF_WIN)
                r[f"N_{l}"] = rho[m].mean()
        if "d" in d:
            E, rho = d["d"]
            r["e_d_occ"] = band_centre(E, rho, E.min(), ef) - ef
        if "p" in d:
            E, rho = d["p"]
            r["e_p_occ"] = band_centre(E, rho, E.min(), ef) - ef
            r["e_p_win"] = band_centre(E, rho, ef + WIN_LO, ef + WIN_HI) - ef
        rows.append(r)

    # ---- consistency check: PDOS integral vs Lowdin charge (first 3 atoms)
    print("\nConsistency check, integral of PDOS up to E_F vs Lowdin charge:")
    chk = [i + 1 for i, e in enumerate(el) if e == "Mg"][:3]     # Mg only: Zn/Ca Lowdin charges include semicore shells outside the PDOS window
    for i in chk:
        for l in ("s", "p"):
            if i in pdos and l in pdos[i]:
                E, rho = pdos[i][l]
                dE = E[1] - E[0]
                integ = rho[E <= ef].sum() * dE
                ql = q.get(i, {}).get(l, np.nan)
                print(f"  atom {i} {l}: integral = {integ:.3f}   Lowdin = {ql:.3f}   ratio = {integ/ql:.2f}")
    E0 = pdos[chk[0]]["s"][0] if chk else next(iter(pdos.values()))["s"][0]
    print(f"  (ratio ~1: PDOS counts both spins; ratio ~0.5: states per spin. Energy step {E0[1]-E0[0]:.3f} eV)")
    top = parse_top_band(f"qe_inputs/{name}.out")
    if top is not None:
        print(f"Highest Kohn-Sham band: {top:.2f} eV = E_F + {top-ef:.2f} eV "
              f"(PDOS has no states above this; broadening adds ~0.15 eV tail)")
        if ef + WIN_HI > top - 0.15:
            print(f"  WARNING: e_p_win upper edge (E_F + {WIN_HI}) is at or beyond the last band -> "
                  f"e_p_win is TRUNCATED. Use e_p_occ and N(E_F), or rerun with more bands (nbnd).")

    # ---- dq relative to the interior of the same patch
    for e in sorted(set(el)):
        ref = [r["q"] for r in rows if r["el"] == e and r["depth"] > INTERIOR_DEPTH and not r["edge"] and not r["bottom"]]
        mean_ref = np.mean(ref) if ref else np.nan   # NaN if this element has no interior atoms in the patch
        # (Zn/Ca in thin interface cells): dq is then undefined - compare q of these atoms across cells instead
        for r in rows:
            if r["el"] == e:
                r["dq"] = r["q"] - mean_ref
                r["dZ"] = r["q"] - ZVAL.get(e, np.nan)     # + = gained electrons (Lowdin; Mg interior ~ 2.00, so meaningful for Mg)

    # ---- grouped summary
    layers = [(f"surface   (<{SURFACE_DEPTH:g} A)", lambda d: d < SURFACE_DEPTH),
              (f"sub-surf  ({SURFACE_DEPTH:g}-{SUBSURFACE_DEPTH:g} A)", lambda d: SURFACE_DEPTH <= d < SUBSURFACE_DEPTH),
              (f"interior  (>{SUBSURFACE_DEPTH:g} A)", lambda d: d >= SUBSURFACE_DEPTH)]
    keys = ["q", "dq", "N_s", "N_p", "e_p_occ", "e_p_win"]
    n_edge = sum(r["edge"] for r in rows)
    n_bot = sum(r["bottom"] for r in rows if not r["edge"])
    print(f"\n{n_edge} of {n} atoms lie within {edge_margin:g} A of the lateral patch edge and are EXCLUDED below, as are {n_bot} more non-edge atoms in the bottom {BOTTOM_MARGIN:g} A (artificial face). All kept in the CSV.")
    print("mean +- std  (n atoms)   [q in e, N in states/eV, e_p in eV rel. E_F]")
    for e in sorted(set(el)):
        for lab, f in layers:
            sel = [r for r in rows if r["el"] == e and f(r["depth"]) and not r["edge"] and not r["bottom"]]
            if not sel:
                continue
            parts = []
            for k in keys:
                v = np.array([r.get(k, np.nan) for r in sel], float)
                parts.append(f"{k}={np.nanmean(v):.3f}+-{np.nanstd(v):.3f}")
            print(f"  {e:2s} {lab:24s} n={len(sel):3d}  " + "  ".join(parts))

    # ---- Mg grouped by solute neighbours: separates the chemical (Zn/Ca) effect from the depth effect
    mg = [r for r in rows if r["el"] == "Mg" and not r["edge"] and not r["bottom"]]
    if mg and any(r["n_ZnCa"] for r in mg):
        print("\nMg atoms grouped by number of Zn/Ca neighbours within %.1f A (all depths, bottom face excluded):" % NN_CUT)
        for k in range(0, 4):
            sel = [r for r in mg if (r["n_ZnCa"] == k if k < 3 else r["n_ZnCa"] >= 3)]
            if len(sel) < 2:
                continue
            f = lambda c: (np.nanmean([r[c] for r in sel]), np.nanstd([r[c] for r in sel]))
            lab = f"{k}" if k < 3 else ">=3"
            print(f"  n_ZnCa={lab:3s} n={len(sel):3d}  q={f('q')[0]:.3f}+-{f('q')[1]:.3f}  e_p_occ={f('e_p_occ')[0]:.3f}+-{f('e_p_occ')[1]:.3f}  N_p={f('N_p')[0]:.3f}")
        qq = np.array([r["q"] for r in mg]); nz = np.array([r["n_ZnCa"] for r in mg], float)
        if nz.std() > 0:
            print(f"  corr(n_ZnCa, q) = {np.corrcoef(nz, qq)[0, 1]:.2f}")
        # surface Mg only, split by neighbours
        sur = [r for r in mg if r["depth"] < SURFACE_DEPTH]
        for lab, sel in (("surface Mg, no Zn/Ca neighbour", [r for r in sur if r["n_ZnCa"] == 0]),
                         ("surface Mg, >=1 Zn/Ca neighbour", [r for r in sur if r["n_ZnCa"] >= 1])):
            if len(sel) >= 2:
                print(f"  {lab:34s} n={len(sel):3d}  q={np.mean([r['q'] for r in sel]):.3f}  e_p_occ={np.mean([r['e_p_occ'] for r in sel]):.3f}")
    zn = [r for r in rows if r["el"] != "Mg" and not r["edge"] and not r["bottom"]]
    for e in sorted({r["el"] for r in zn}):
        sel = [r for r in zn if r["el"] == e]
        ed = [r["e_d_occ"] for r in sel if "e_d_occ" in r]
        print(f"  {e}: n={len(sel)}  mean q={np.mean([r['q'] for r in sel]):.3f}  dZ={np.mean([r['dZ'] for r in sel]):+.3f} e"
              + (f"  d-band centre (occupied)={np.mean(ed):.2f} eV" if ed else "  (no d PDOS found)"))

    # ---- diagnostic only: the artificial bottom face (excluded above). In a symmetric slab it should
    #      match the top surface, which is a built-in check of the whole pipeline.
    for e in sorted(set(el)):
        sel = [r for r in rows if r["el"] == e and r["bottom"] and not r["edge"]]
        if sel:
            v = lambda k: np.nanmean([r.get(k, np.nan) for r in sel])
            print(f"  {e:2s} bottom face (diagnostic only) n={len(sel):3d}  q={v('q'):.3f}  e_p_occ={v('e_p_occ'):.3f}")

    # ---- csv
    cols = ["atom", "el", "x", "y", "z", "depth", "edge", "bottom", "n_ZnCa", "q", "q_s", "q_p", "dq", "dZ", "N_s", "N_p", "e_p_occ", "e_p_win", "e_d_occ"]
    with open(f"analysis_{name}.csv", "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(str(r[c]) if c in ("atom", "el", "edge", "bottom", "n_ZnCa") else
                             ("%.5f" % r[c] if c in r else "")
                             for c in cols) + "\n")
    print(f"\nwrote analysis_{name}.csv")


if __name__ == "__main__":
    # usage: 06_analyze_pdos.py NAME [geometry.xyz]   (2nd arg for k-point-tagged runs, e.g.
    #        python3 06_analyze_pdos.py crystal_ideal_k221 slab_crystal_ideal.xyz)
    main(sys.argv[1] if len(sys.argv) > 1 else "Mg60Zn35Ca5_crystal_core",
         sys.argv[2] if len(sys.argv) > 2 else None)
