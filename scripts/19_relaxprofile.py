"""
Does Zn really relax outward, or does the 1.5 A window just shrink?  (no DFT, seconds)

18_segtest.py showed the whole outermost-layer Zn excess is made by the MEAM relaxation, not the cut.
Two mechanisms give that signature, and a fixed DEPTH window cannot tell them apart:
  (a) differential outward relaxation - Zn moves out relative to Mg (real in the model)
  (b) roughening - the window is set by the single highest atom, so a rougher surface samples fewer
      atoms and whichever element sticks out dominates them (pure sampling, no differential motion)
Mean exposed count falls 18.5 -> 13.4 on relaxation, so (b) is live.

Three tests, glass patches only (iface_*/crystal_* excluded: a Mg-only nanocrystal makes the
whole-patch baseline meaningless, e.g. iface_N_D26 reads +94 points purely from dilution):
  1. FIXED COUNT  topmost N atoms in both, same N -> sampling identical by construction,
                  so only real displacement can move the composition. Run at several N.
  2. ROUGHNESS    spread of the top 20 atoms, and how many fall in 1.5 A, cut vs relaxed.
  3. HEIGHT       per-element mean height above the slab median, cut vs relaxed. If atom order
                  survives the relax, this is done per atom and is the direct measure.

  python3 19_relaxprofile.py
  KEEP_IFACE=1 python3 19_relaxprofile.py     # include the inhomogeneous patches
"""
import glob, os, sys
import numpy as np

ELS = ("Ca", "Mg", "Zn")
TYPE = {1: "Mg", 2: "Zn", 3: "Ca"}
COUNTS = [int(x) for x in os.environ.get("COUNTS", "10 15 20 25 30").split()]
KEEP = os.environ.get("KEEP_IFACE", "") == "1"
rng = np.random.default_rng(0)


def read_data(fn):
    L = open(fn).read().splitlines()
    n = int(next(l for l in L if l.strip().endswith("atoms")).split()[0])
    i = next(k for k, l in enumerate(L) if l.strip().startswith("Atoms")) + 2
    rec = []
    for l in L[i:i + n]:
        p = l.split()
        rec.append((int(p[0]), TYPE.get(int(p[1]), "?"), float(p[4])))
    rec.sort()                                   # by atom id
    return np.array([r[1] for r in rec]), np.array([r[2] for r in rec])


def read_xyz(fn):
    L = open(fn).read().splitlines()
    n = int(L[0].split()[0])
    at = [l.split() for l in L[2:2 + n]]
    return np.array([a[0] for a in at]), np.array([float(a[3]) for a in at])


def topn(el, z, n):
    k = np.argsort(-z)[:n]
    return el[k]


def perm(d):
    d = np.asarray(d, float)
    d = d[~np.isnan(d)]
    if len(d) < 3:
        return np.nan, np.nan, np.nan, 0
    m = d.mean()
    c = sum(abs((d * rng.choice([-1, 1], len(d))).mean()) >= abs(m) - 1e-12 for _ in range(20000))
    return m, d.std(ddof=1) / np.sqrt(len(d)), (c + 1) / 20001, len(d)


names = sorted(os.path.basename(f)[5:-4] for f in glob.glob("slab_*.xyz")
               if os.path.exists("patches/" + os.path.basename(f)[5:-4] + ".data"))
if not KEEP:
    names = [n for n in names if not n.startswith(("iface_", "crystal_"))]
rows = []
for nm in names:
    try:
        eu, zu = read_data(f"patches/{nm}.data")
        er, zr = read_xyz(f"slab_{nm}.xyz")
    except Exception as ex:
        print(f"skip {nm}: {ex}"); continue
    rows.append((nm, eu, zu, er, zr))
if not rows:
    print("no patch pairs found"); sys.exit()
print(f"{len(rows)} glass patches" + ("" if not KEEP else " (iface/crystal included)") + "\n")

# ---- 1. fixed count -------------------------------------------------------
print("1. FIXED COUNT  Zn% among the topmost N atoms, same N in both (sampling held constant)")
print(f"   {'N':>4s} {'cut Zn%':>8s} {'relaxed Zn%':>12s} {'difference':>12s} {'se':>5s} {'p':>7s}")
for N in COUNTS:
    d, cu, cr = [], [], []
    for nm, eu, zu, er, zr in rows:
        if min(len(eu), len(er)) < N:
            continue
        a = 100 * (topn(eu, zu, N) == "Zn").mean()
        b = 100 * (topn(er, zr, N) == "Zn").mean()
        cu.append(a); cr.append(b); d.append(b - a)
    m, se, p, k = perm(d)
    print(f"   {N:4d} {np.mean(cu):8.1f} {np.mean(cr):12.1f} {m:+12.1f} {se:5.1f} {p:7.4f}   ({k} patches)")

# ---- 2. roughness ---------------------------------------------------------
print("\n2. ROUGHNESS  of the surface, cut vs relaxed")
sp, nw = [], []
for nm, eu, zu, er, zr in rows:
    a = np.sort(zu)[-20:]; b = np.sort(zr)[-20:]
    sp.append((b.max() - b.min()) - (a.max() - a.min()))
    na = ((zu.max() - zu < 1.5) & (zu - zu.min() >= 2.0)).sum()
    nb = ((zr.max() - zr < 1.5) & (zr - zr.min() >= 2.0)).sum()
    nw.append(nb - na)
for lab, d, u in (("spread of the top 20 atoms", sp, "A"), ("atoms inside the 1.5 A window", nw, "atoms")):
    m, se, p, k = perm(d)
    print(f"   {lab:32s} {m:+6.2f} +- {se:4.2f} {u:6s} p = {p:.4f}")

# ---- 3. per-element height ------------------------------------------------
matched = all(len(eu) == len(er) and (eu == er).all() for _, eu, zu, er, zr in rows)
print("\n3. HEIGHT above the slab median, per element, cut vs relaxed"
      + ("   [atom order preserved -> per-atom]" if matched else "   [order not preserved -> distributional]"))
print(f"   {'element':>8s} {'cut':>7s} {'relaxed':>8s} {'moved':>7s} {'se':>5s} {'p':>7s}")
mv = {}
for e in ELS:
    d, hu_, hr_ = [], [], []
    for nm, eu, zu, er, zr in rows:
        hu = zu - np.median(zu); hr = zr - np.median(zr)
        if matched:
            m_ = (eu == e) & (hu > np.percentile(hu, 90))   # near-surface atoms of this element
            if m_.sum() < 2:
                continue
            a, b = hu[m_].mean(), hr[m_].mean()
        else:
            mu = (eu == e) & (hu > np.percentile(hu, 90))
            mr = (er == e) & (hr > np.percentile(hr, 90))
            if mu.sum() < 2 or mr.sum() < 2:
                continue
            a, b = hu[mu].mean(), hr[mr].mean()
        hu_.append(a); hr_.append(b); d.append(b - a)
    m, se, p, k = perm(d)
    mv[e] = (m, se)
    print(f"   {e:>8s} {np.mean(hu_):7.2f} {np.mean(hr_):8.2f} {m:+7.2f} {se:5.2f} {p:7.4f}   ({k} patches)")
if "Zn" in mv and "Mg" in mv and not np.isnan(mv["Zn"][0]):
    dz = mv["Zn"][0] - mv["Mg"][0]
    print(f"\n   Zn relative to Mg: {dz:+.2f} A"
          + ("  -> Zn relaxes outward" if dz > 0.02 else
             "  -> no differential outward motion of Zn"))

print("""
Reading it:
  FIXED COUNT significant and positive   -> real differential outward relaxation of Zn.
                                            Report as relaxation, not segregation, and flag that this
                                            potential underbinds Zn and so would exaggerate it.
  FIXED COUNT ~0 while the window shrank -> the 1.5 A excess is a sampling artifact of a roughened
                                            surface. Drop the claim; the exposed fractions the binding
                                            model uses stay valid, since the DFT sites sit on that
                                            same relaxed surface either way.""")
