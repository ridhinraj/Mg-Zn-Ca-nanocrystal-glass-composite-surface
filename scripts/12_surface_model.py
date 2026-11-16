"""
Surface-composition model of water reactivity (no DFT; reads CSVs only).

Inputs : analysis_<cell>.csv (every cell, from 06_analyze_pdos.py)   ads/ads_table.csv (from 11_adsorption.py analyze)
Model A: E_pred(cell) = sum_el f_el(cell) * <E_matched(el)>       f_el = fraction of exposed top atoms (depth < DEPTH_MAX)
Model B: as A, but each exposed atom gets <E(el)> + slope_el * (dZ_atom - <dZ(el)>)   (per-element dZ slope, only if |rho|>0.4)
Checks : leave-one-site-out error of A, B and a region-only model (does the element-resolved model beat region means?)
Output : surface_model.csv (per cell), printed table with bootstrap 95% intervals (resampling the DFT sites per element).
Physisorbed fraction: fraction of exposed atoms that are Zn x P(physisorbed | Zn) from the DFT sites (CN(O)==0).
"""
import csv, glob, os, re
import numpy as np

DEPTH_MAX = float(os.environ.get("DEPTH_MAX", "1.5"))
ZVAL = {"Mg": 2.0, "Zn": 20.0, "Ca": 10.0}
ELS = ("Ca", "Mg", "Zn")
rng = np.random.default_rng(3)

T = [r for r in csv.DictReader(open("ads/ads_table.csv")) if r["mol"] == "H2O" and r["cut"] == "std" and r["Em"] not in ("nan", "")]
T = [r for r in T if np.isfinite(float(r["Em"])) and r["dZ"] != ""]
E = {el: np.array([float(r["Em"]) for r in T if r["el"] == el]) for el in ELS}
DZ = {el: np.array([float(r["dZ"]) for r in T if r["el"] == el]) for el in ELS}
PHYS = {el: np.mean([int(r["cn"]) == 0 for r in T if r["el"] == el]) for el in ELS}
print("DFT sites:", {el: len(E[el]) for el in ELS}, " mean E:", {el: round(E[el].mean(), 3) for el in ELS})


def slope(el):
    x, y = DZ[el], E[el]
    if len(x) < 5 or x.std() == 0:
        return 0.0
    rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
    return float(np.polyfit(x, y, 1)[0]) if abs(np.corrcoef(rx, ry)[0, 1]) > 0.4 else 0.0


SL = {el: slope(el) for el in ELS}
print("dZ slopes used (eV per e):", {el: round(v, 3) for el, v in SL.items()})

cells = sorted(os.path.basename(f)[9:-4] for f in glob.glob("analysis_*.csv"))
out = []
# ---- cells known only from geometry (no PDOS needed): EXTRA_IN="qe_inputs/A_wf.in qe_inputs/B_wf.in" -> exposed fractions
def geom_counts(fn):
    t = open(fn).read().split("ATOMIC_POSITIONS angstrom", 1)[1].split("K_POINTS", 1)[0]
    at = [l.split() for l in t.strip().splitlines() if len(l.split()) >= 4]
    el = np.array([a[0] for a in at]); z = np.array([float(a[3]) for a in at])
    m = (z.max() - z < DEPTH_MAX) & (z - z.min() >= 2.0)
    return {e: int((el[m] == e).sum()) for e in ELS}, int(m.sum())


def geom_fractions(fn):
    c, n = geom_counts(fn)
    return ({e: c[e] / n for e in ELS} if n else {e: 0.0 for e in ELS}), n


def predict(f):
    return sum(f[e] * E[e].mean() for e in ELS)
series = {}
for fn in os.environ.get("EXTRA_IN", "").split():
    cnt, nn = geom_counts(fn)
    if not nn:
        print(f"skipping {fn} (no exposed atoms)"); continue
    f = {e: cnt[e] / nn for e in ELS}
    A = predict(f)
    bs = [sum(f[e] * rng.choice(E[e], len(E[e])).mean() for e in ELS) for _ in range(4000)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    name = os.path.basename(fn)[:-3]
    out.append(dict(cell=name, n=nn, fCa=f["Ca"], fMg=f["Mg"], fZn=f["Zn"], E_A=A, lo=lo, hi=hi, E_B=float("nan"),
                    phys=sum(f[e] * PHYS[e] for e in ELS)))
    m = re.match(r"(c\d+_[a-z]+\d+)_", name.replace("_wf", ""))
    if m:
        series.setdefault(m.group(1), []).append((cnt, nn, A))
for c in cells:
    rd = csv.DictReader(open(f"analysis_{c}.csv"))
    if not rd.fieldnames or not {"bottom", "depth", "el"} <= set(rd.fieldnames):
        print(f"skipping analysis_{c}.csv (columns: {rd.fieldnames})"); continue
    rows = [r for r in rd if r["bottom"] and int(float(r["bottom"])) == 0 and r["depth"] and float(r["depth"]) < DEPTH_MAX]
    if not rows:
        print(f"skipping analysis_{c}.csv (no exposed atoms)"); continue
    el = np.array([r["el"] for r in rows])
    dz = np.array([float(r["q"]) - ZVAL[r["el"]] if r.get("q", "") not in ("", "nan") else np.nan for r in rows])
    f = {e: float((el == e).mean()) for e in ELS}
    A = sum(f[e] * E[e].mean() for e in ELS)
    B = 0.0
    for e in ELS:
        m = el == e
        if m.any():
            d = dz[m]; d = d[np.isfinite(d)]
            B += f[e] * (E[e].mean() + (SL[e] * (d.mean() - DZ[e].mean()) if len(d) else 0.0))
    bs = []
    for _ in range(4000):
        bs.append(sum(f[e] * rng.choice(E[e], len(E[e])).mean() for e in ELS))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    phys = sum(f[e] * PHYS[e] for e in ELS)
    out.append(dict(cell=c, n=len(rows), fCa=f["Ca"], fMg=f["Mg"], fZn=f["Zn"], E_A=A, lo=lo, hi=hi, E_B=B, phys=phys))

print(f"\n{'cell':16s} {'n':>3s} {'fCa':>5s} {'fMg':>5s} {'fZn':>5s}  {'E_A/eV':>7s} [95% CI]          {'E_B/eV':>7s} f_physisorbed")
for r in out:
    print(f"{r['cell']:16s} {r['n']:3d} {r['fCa']:5.2f} {r['fMg']:5.2f} {r['fZn']:5.2f}  {r['E_A']:+7.3f} [{r['lo']:+.3f},{r['hi']:+.3f}]  {r['E_B']:+7.3f} {r['phys']:.2f}")
with open("surface_model.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)

# ---- leave-one-site-out validation on the DFT sites
reg = lambda c: "crystal" if c.startswith("crystal") else ("interphase" if c.startswith("iface") else "glass")
y = np.array([float(r["Em"]) for r in T]); el_ = np.array([r["el"] for r in T]); rg = np.array([reg(r["cell"]) for r in T])
dzs = np.array([float(r["dZ"]) for r in T])
def loo(kind):
    err = []
    for i in range(len(y)):
        m = np.arange(len(y)) != i
        if kind == "grand":
            p = y[m].mean()
        elif kind == "region":
            p = y[m & (rg == rg[i])].mean() if (m & (rg == rg[i])).any() else y[m].mean()
        elif kind == "element":
            p = y[m & (el_ == el_[i])].mean()
        else:  # element + dZ slope
            mm = m & (el_ == el_[i])
            x = dzs[mm]
            p = y[mm].mean() if x.std() == 0 or mm.sum() < 5 else np.polyval(np.polyfit(x, y[mm], 1), dzs[i])
        err.append(y[i] - p)
    err = np.array(err); return np.sqrt((err ** 2).mean())
print("\nLeave-one-site-out RMSE of E_ads (eV), lower = better:")
for k in ("grand", "region", "element", "element+dZ"):
    print(f"  {k:11s} {loo(k):.3f}")


# ---- (1) crystal check: predicted pure-Mg value vs DFT site-matched H2O on the crystal slab (out of sample for the amorphous fit)
C = [r for r in csv.DictReader(open("ads/ads_table.csv")) if r["cell"].startswith("crystal") and r["mol"] == "H2O" and r["cut"] == "std"]
print("\nCrystal check (Mg(0001), H2O, site-matched). Model A predicts the Mg-site mean from the amorphous cells:")
print(f"  predicted pure-Mg: {E['Mg'].mean():+.3f} eV (spread of the amorphous Mg sites +-{E['Mg'].std(ddof=1):.3f})")
for r in C:
    print(f"  {r['cell']:20s} {r['site']:5s} E_matched = {float(r['Em']):+.3f}  E_unmatched = {float(r['Eads']):+.3f}  CN(O)={r['cn']}")

# ---- (2) leave-one-CELL-out: predict each DFT cell's own element-weighted mean from element means fitted on the other cells
print("\nLeave-one-cell-out: cell-level measured vs predicted (cells with >= 3 DFT sites)")
cl = np.array([r["cell"] for r in T]); pr, ms, lab = [], [], []
for c in sorted(set(cl)):
    m = cl == c
    if m.sum() < 3 or c.startswith("crystal"):
        continue
    f = {e: float((el_[m] == e).mean()) for e in ELS}              # element mix of the DFT sites of this cell
    ok = all(((el_ == e) & ~m).sum() > 0 for e in ELS if f[e] > 0)
    if not ok:
        continue
    pred = sum(f[e] * y[(el_ == e) & ~m].mean() for e in ELS if f[e] > 0)
    pr.append(pred); ms.append(y[m].mean()); lab.append(c)
    print(f"  {c:14s} n={m.sum():2d}  measured {y[m].mean():+.3f}  predicted {pred:+.3f}")
if len(pr) > 2:
    pr, ms = np.array(pr), np.array(ms)
    print(f"  RMSE {np.sqrt(((pr - ms) ** 2).mean()):.3f} eV over {len(pr)} cells (cell means of stratified samples: indicative only)")


# ---- composition series: patches of one composition pooled, so the exposed fractions rest on ~240 atoms
if series:
    print("\n" + "=" * 78)
    print("COMPOSITION SERIES, patches pooled per composition")
    print(f"{'composition':14s} {'patches':>7s} {'exposed':>7s} {'fCa':>5s} {'fMg':>5s} {'fZn':>5s}  {'E_pred/eV':>9s}  "
          f"{'patch sd':>8s} {'s.e.':>6s}  [95% CI incl. DFT sites]")
    rows = []
    for lab in sorted(series):
        ps = series[lab]
        tot = {e: sum(c[e] for c, _, _ in ps) for e in ELS}; N = sum(n for _, n, _ in ps)
        f = {e: tot[e] / N for e in ELS}
        Ep = predict(f)
        per = np.array([a for _, _, a in ps])
        sd = per.std(ddof=1) if len(per) > 1 else float("nan")
        se = sd / np.sqrt(len(per)) if len(per) > 1 else float("nan")
        bs = []
        for _ in range(4000):
            k = rng.integers(0, len(ps), len(ps))
            t2 = {e: sum(ps[j][0][e] for j in k) for e in ELS}; n2 = sum(ps[j][1] for j in k)
            f2 = {e: t2[e] / n2 for e in ELS}
            bs.append(sum(f2[e] * rng.choice(E[e], len(E[e])).mean() for e in ELS))
        lo, hi = np.percentile(bs, [2.5, 97.5])
        rows.append((lab, f, Ep, sd, se, lo, hi, len(ps), N))
        print(f"{lab:14s} {len(ps):7d} {N:7d} {f['Ca']:5.2f} {f['Mg']:5.2f} {f['Zn']:5.2f}  {Ep:+9.3f}  "
              f"{sd:8.3f} {se:6.3f}  [{lo:+.3f},{hi:+.3f}]")
    with open("series_model.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["composition", "patches", "exposed", "fCa", "fMg", "fZn", "E_pred", "patch_sd", "se", "lo", "hi"])
        for lab, f, Ep, sd, se, lo, hi, np_, N in rows:
            w.writerow([lab, np_, N, f["Ca"], f["Mg"], f["Zn"], Ep, sd, se, lo, hi])
    if len(rows) > 2:
        sp = max(r[2] for r in rows) - min(r[2] for r in rows)
        mse = np.nanmean([r[4] for r in rows])
        print(f"\n  spread across compositions {sp:.3f} eV   mean standard error {mse:.3f} eV   ratio {sp / mse:.1f}")
        print("  (a ratio well above 2 means the composition trend is resolved; near 1 means more patches are needed)")
