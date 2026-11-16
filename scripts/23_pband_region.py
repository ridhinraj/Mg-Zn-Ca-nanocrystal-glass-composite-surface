"""
Is the p-band/neighbour slope really different between interphase and glass matrix?  (23_pband_region.py)

Fig. 2 shows every interphase cell with a POSITIVE slope and every glass-matrix cell NEGATIVE, which
would be a region difference in a paper whose headline is that the two regions are indistinguishable.
Dropping the zero-neighbour crystal Mg did not remove it. Two things could still produce it without any
real region effect:
  * the populations barely overlap - interphase Mg sit at 1-3 Zn/Ca neighbours, glass Mg at 2-6 - so two
    straight lines are being fitted to different windows of one possibly curved relationship
  * the slopes may not be individually significant at all: the scatter is about +-0.25 eV

Three tests:
  1. SLOPES with a permutation p per cell, on all atoms, on n>=1, and on the common window.
  2. MATCHED COUNT - the decisive one. At the SAME neighbour count, do interphase and glass-matrix Mg
     have the same p-band centre? No window, no line fitting, no extrapolation.
  3. POOLED SHAPE - mean p-band against neighbour count for each region, to see whether one curve
     through both regions explains the data.

Cells within a region are NOT independent (three cut depths at one pole share glass atoms), so the
region-level test is also run on pole/cell groups, which is the honest unit.

  python3 23_pband_region.py
"""
import csv, glob, os, re
import numpy as np

rng = np.random.default_rng(0)
MINN = int(os.environ.get("MINN", "8"))


def region(cell):
    if cell.startswith("iface_"):
        return "interphase"
    if cell.startswith("matrix_"):
        return "glass matrix"
    return "crystal" if cell.startswith("crystal") else "other"


def group(cell):
    """independent realisation: the three depths at one pole share glass atoms"""
    m = re.match(r"iface_([NS])_D\d+", cell)
    return f"iface_{m.group(1)}" if m else cell


def load():
    out = {}
    for fn in sorted(glob.glob("analysis_*.csv")):
        cell = os.path.basename(fn)[9:-4]
        if region(cell) not in ("interphase", "glass matrix"):
            continue
        rr = list(csv.DictReader(open(fn)))
        if not rr or "n_ZnCa" not in rr[0] or "e_p_win" not in rr[0]:
            continue
        p = []
        for r in rr:
            if r.get("el") != "Mg":
                continue
            try:
                p.append((float(r["n_ZnCa"]), float(r["e_p_win"])))
            except (TypeError, ValueError):
                pass
        p = np.array([q for q in p if np.isfinite(q).all()])
        if len(p) >= MINN:
            out[cell] = p
    return out


def slope_p(x, y, n=5000):
    if len(set(x)) < 3 or len(x) < 6:
        return float("nan"), float("nan"), float("nan")
    k, b = np.polyfit(x, y, 1)
    r2 = 1 - np.sum((y - (k * x + b)) ** 2) / np.sum((y - y.mean()) ** 2)
    c = sum(abs(np.polyfit(x, rng.permutation(y), 1)[0]) >= abs(k) - 1e-15 for _ in range(n))
    return k, r2, (c + 1) / (n + 1)


D = load()
if not D:
    raise SystemExit("no analysis_*.csv with Mg n_ZnCa and e_p_win")
allx = np.concatenate([p[:, 0] for p in D.values()])
byreg = {}
for c, p in D.items():
    byreg.setdefault(region(c), []).append((c, p))

# common window: neighbour counts where BOTH regions have a decent number of atoms
cnt = {r: {} for r in byreg}
for r, lst in byreg.items():
    for _, p in lst:
        for v in p[:, 0]:
            cnt[r][v] = cnt[r].get(v, 0) + 1
common = sorted(v for v in set(cnt["interphase"]) & set(cnt["glass matrix"])
                if min(cnt["interphase"][v], cnt["glass matrix"][v]) >= 15)
print(f"{len(D)} cells, {len(allx)} Mg atoms.  Common neighbour counts (>=15 atoms in each region): "
      + (", ".join(f"{int(v)}" for v in common) if common else "none"))

print("\n1. SLOPES  (eV per neighbour; p from permuting the p-band values within the cell)")
print(f"   {'cell':16s} {'region':13s} {'all atoms':>22s} {'n>=1':>22s} {'common window':>22s}")
agg = {}
for r in ("interphase", "glass matrix"):
    for c, p in sorted(byreg.get(r, [])):
        row = []
        for sel in (p, p[p[:, 0] >= 1], p[np.isin(p[:, 0], common)] if common else p[:0]):
            k, r2, pv = slope_p(sel[:, 0], sel[:, 1]) if len(sel) else (np.nan,) * 3
            row.append((k, r2, pv))
        agg.setdefault(r, []).append((c, row))
        print(f"   {c:16s} {r:13s} " + " ".join(
            f"{k:+.3f} R2={r2:4.2f} p={pv:5.3f}" if np.isfinite(k) else f"{'-':>22s}" for k, r2, pv in row))
print("\n   region means")
for r, v in agg.items():
    for j, lab in enumerate(("all atoms", "n>=1", "common window")):
        ks = [row[j][0] for _, row in v if np.isfinite(row[j][0])]
        sig = sum(1 for _, row in v if np.isfinite(row[j][2]) and row[j][2] < 0.05)
        if ks:
            print(f"     {r:13s} {lab:15s} {np.mean(ks):+.3f}  ({sig}/{len(ks)} cells with p < 0.05)")

print("\n2. MATCHED NEIGHBOUR COUNT  -  p-band centre of interphase minus glass matrix, at the same n")
print(f"   {'n':>3s} {'interphase':>22s} {'glass matrix':>22s} {'difference':>12s} {'p':>7s}")
diffs = []
for v in common:
    a = np.concatenate([p[p[:, 0] == v, 1] for _, p in byreg["interphase"]])
    b = np.concatenate([p[p[:, 0] == v, 1] for _, p in byreg["glass matrix"]])
    obs = a.mean() - b.mean()
    pool = np.concatenate([a, b]); k = len(a)
    c = sum(abs((q := rng.permutation(pool))[:k].mean() - q[k:].mean()) >= abs(obs) - 1e-15 for _ in range(5000))
    diffs.append((v, obs))
    print(f"   {int(v):3d} {a.mean():+9.3f} +- {a.std(ddof=1) / np.sqrt(len(a)):.3f} (n={len(a):4d}) "
          f"{b.mean():+9.3f} +- {b.std(ddof=1) / np.sqrt(len(b)):.3f} (n={len(b):4d}) "
          f"{obs:+12.3f} {(c + 1) / 5001:7.4f}")
if diffs:
    d = np.array([x for _, x in diffs])
    print(f"\n   mean difference across matched counts: {d.mean():+.3f} eV"
          + ("   -> the regions DIFFER at matched neighbour count: a real region effect"
             if abs(d.mean()) > 0.05 else
             "   -> at matched neighbour count the regions AGREE: the slope flip is a range artifact"))

print("\n3. POOLED SHAPE  -  mean p-band at each neighbour count, per region")
print(f"   {'n':>3s} " + "".join(f"{r:>18s}" for r in ("interphase", "glass matrix")))
for v in sorted(set(allx)):
    cells = []
    for r in ("interphase", "glass matrix"):
        y = [p[p[:, 0] == v, 1] for _, p in byreg.get(r, [])]
        y = np.concatenate(y) if y else np.array([])
        cells.append(f"{y.mean():+8.3f} (n={len(y):4d})" if len(y) >= 5 else f"{'-':>18s}")
    print(f"   {int(v):3d} " + "".join(f"{c:>18s}" for c in cells))

print("\n4. INDEPENDENT UNITS  (three depths at one pole share glass atoms, so cells are not independent)")
gslopes = {}
for r, v in agg.items():
    for c, row in v:
        if np.isfinite(row[1][0]):
            gslopes.setdefault((r, group(c)), []).append(row[1][0])
units = {k: np.mean(x) for k, x in gslopes.items()}
A = [v for (r, _), v in units.items() if r == "interphase"]
B = [v for (r, _), v in units.items() if r == "glass matrix"]
print(f"   interphase units   {len(A)}: " + ", ".join(f"{v:+.3f}" for v in A))
print(f"   glass matrix units {len(B)}: " + ", ".join(f"{v:+.3f}" for v in B))
if A and B:
    obs = abs(np.mean(A) - np.mean(B)); pool = np.array(A + B); k = len(A)
    c = sum(abs((q := rng.permutation(pool))[:k].mean() - q[k:].mean()) >= obs - 1e-15 for _ in range(20000))
    print(f"   difference {np.mean(A) - np.mean(B):+.3f} eV per neighbour, permutation p = {(c + 1) / 20001:.4f} "
          f"on {len(A)} vs {len(B)} independent units")
    print("   (with so few units this test has little power: read it alongside test 2, which uses atoms)")
print("""
How to read it: test 2 decides. If interphase and glass-matrix Mg agree at matched neighbour count, the
opposite slopes come from fitting lines to two non-overlapping windows of one curve, and Fig. 2 should
show the pooled relationship instead of per-cell lines. If they differ at matched count, there IS a
region effect in the electronic structure, and the Discussion has to say so and reconcile it with the
work-function and binding nulls - which are direct measurements and are not affected by this.""")
