"""
Is the outermost-layer Zn excess made by the CUT or by the RELAXATION?  (no DFT, seconds)

For every patch it compares the SAME depth window on two structures:
  patches/<name>.data   the cut, before MEAM relaxation
  slab_<name>.xyz       after relaxation (what 12_surface_model.py actually reads)
Exposed set in both: z within DEPTH_MAX (1.5 A) of the highest atom and at least 2.0 A above the lowest.

  unrelaxed excess  = cut composition is already Zn-rich -> the cut or the exposed-set definition does it
  relaxed - unrelax = relaxation moves Zn outward      -> differential outward relaxation, a property of the potential
Zn cannot diffuse in the 50 ps relaxation, so a large second term is NOT Gibbsian segregation and would be
sensitive to the potential's known underbinding of Zn.

  python3 18_segtest.py                 # every patch found
  PATCHES="c04_zn25_60_0 ..." python3 18_segtest.py
"""
import glob, os, sys
import numpy as np

DEPTH_MAX = float(os.environ.get("DEPTH_MAX", "1.5"))
ELS = ("Ca", "Mg", "Zn")
TYPE = {1: "Mg", 2: "Zn", 3: "Ca"}
rng = np.random.default_rng(0)


def read_data(fn):
    L = open(fn).read().splitlines()
    n = int(next(l for l in L if l.strip().endswith("atoms")).split()[0])
    i = next(k for k, l in enumerate(L) if l.strip().startswith("Atoms")) + 2
    el, z = [], []
    for l in L[i:i + n]:
        p = l.split()
        el.append(TYPE.get(int(p[1]), "?")); z.append(float(p[4]))
    return np.array(el), np.array(z)


def read_xyz(fn):
    L = open(fn).read().splitlines()
    n = int(L[0].split()[0])
    at = [l.split() for l in L[2:2 + n]]
    return np.array([a[0] for a in at]), np.array([float(a[3]) for a in at])


def top(el, z):
    m = (z.max() - z < DEPTH_MAX) & (z - z.min() >= 2.0)
    return ({e: int((el[m] == e).sum()) for e in ELS}, int(m.sum()),
            {e: float((el == e).mean()) for e in ELS})


names = os.environ.get("PATCHES", "").split() or sorted(
    os.path.basename(f)[5:-4] for f in glob.glob("slab_*.xyz")
    if os.path.exists("patches/" + os.path.basename(f)[5:-4] + ".data"))
rows = []
for nm in names:
    try:
        eu, zu = read_data(f"patches/{nm}.data")
        er, zr = read_xyz(f"slab_{nm}.xyz")
    except Exception as ex:
        print(f"skip {nm}: {ex}"); continue
    cu, nu, wu = top(eu, zu)
    cr, nr, wr = top(er, zr)
    rows.append((nm, cu, nu, wu, cr, nr, wr))
if not rows:
    print("no patch pairs found (need patches/<name>.data and slab_<name>.xyz)"); sys.exit()

print(f"{len(rows)} patches, outermost {DEPTH_MAX} A\n")
print(f"{'patch':24s} {'whole Zn':>8s} | {'cut: n':>6s} {'Zn%':>5s} {'excess':>7s} | {'relaxed: n':>10s} {'Zn%':>5s} {'excess':>7s} | {'relax moves':>11s}")
A, B = [], []
for nm, cu, nu, wu, cr, nr, wr in rows:
    zu_ = 100 * cu["Zn"] / nu if nu else np.nan
    zr_ = 100 * cr["Zn"] / nr if nr else np.nan
    w = 100 * wr["Zn"]
    A.append(zu_ - w); B.append(zr_ - zu_)
    print(f"{nm:24s} {w:7.1f}% | {nu:6d} {zu_:5.1f} {zu_ - w:+7.1f} | {nr:10d} {zr_:5.1f} {zr_ - w:+7.1f} | {zr_ - zu_:+11.1f}")
A, B = np.array(A), np.array(B)


def perm(d):
    m = d.mean()
    c = sum(abs((d * rng.choice([-1, 1], len(d))).mean()) >= abs(m) - 1e-12 for _ in range(20000))
    return m, d.std(ddof=1) / np.sqrt(len(d)), (c + 1) / 20001


for lab, d in (("CUT alone (unrelaxed exposed - whole patch)", A),
               ("RELAXATION alone (relaxed - unrelaxed, same window)", B),
               ("TOTAL (relaxed exposed - whole patch)", A + B)):
    m, se, p = perm(d)
    print(f"\n{lab}\n  Zn {m:+5.1f} +- {se:4.1f} points   sign-permutation p = {p:.4f}")
print("\nIf CUT is ~0 and RELAXATION carries the whole effect, this is outward relaxation of weakly bound Zn,")
print("not compositional segregation: report it as such and flag the potential's Zn underbinding.")
