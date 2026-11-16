"""
Does local coordination explain water binding beyond the element?  (no DFT; reads ads/ads_table.csv, ads/sites_*.csv, qe_inputs/<cell>_wf.in)

Per adsorption site (atom under the water): total metal coordination CN within 3.4 A (also 3.2 / 3.6 as a cutoff check), depth below the
top atom, Zn/Ca neighbours, and the clean-slab excess charge dZ. Then:
  * rank correlations of E_matched with CN (pooled and within each element, permutation p-values)
  * leave-one-site-out RMSE of  element | element+CN | element+depth | element+CN+dZ | element x CN (own slope per element)
  * transfer test: fit on the amorphous sites, predict the Mg(0001) atop site (DFT -0.381 eV; element-only prediction was -0.569 eV)
Run via submit_coordination.sh (testq, 1 core).
"""
import csv, glob, os
import numpy as np

ZVAL = {"Mg": 2.0, "Zn": 20.0, "Ca": 10.0}
ELS = ("Ca", "Mg", "Zn")
rng = np.random.default_rng(5)


def load_atoms(fn):
    t = open(fn).read()
    cell = [list(map(float, l.split())) for l in t.split("CELL_PARAMETERS angstrom", 1)[1].split("ATOMIC_POSITIONS", 1)[0].strip().splitlines()]
    pos = t.split("ATOMIC_POSITIONS angstrom", 1)[1].split("K_POINTS", 1)[0]
    at = [(l.split()[0], float(l.split()[1]), float(l.split()[2]), float(l.split()[3])) for l in pos.strip().splitlines() if len(l.split()) >= 4]
    return cell[0][0], cell[1][1], at


def final_atoms(fn):
    if not os.path.exists(fn):
        return []
    t = open(fn, errors="ignore").read()
    out = []
    for l in t.split("ATOMIC_POSITIONS")[-1].splitlines()[1:]:
        q = l.split()
        if len(q) < 4:
            break
        try:
            out.append((q[0], float(q[1]), float(q[2]), float(q[3])))
        except ValueError:
            break
    return out


def descriptors(cell, atom1, site):
    Lx, Ly, at = load_atoms(f"qe_inputs/{cell}_wf.in")
    P = np.array([[a[1], a[2], a[3]] for a in at]); el = np.array([a[0] for a in at])
    i = atom1 - 1
    d = P - P[i]
    d[:, 0] -= Lx * np.round(d[:, 0] / Lx); d[:, 1] -= Ly * np.round(d[:, 1] / Ly)
    r = np.sqrt((d ** 2).sum(1)); r[i] = 1e9
    out = {f"cn{int(c * 10)}": int((r < c).sum()) for c in (3.2, 3.4, 3.6)}
    out["depth"] = float(P[:, 2].max() - P[i, 2])
    out["nZnCa"] = int(((r < 3.4) & np.isin(el, ("Zn", "Ca"))).sum())
    # PREDICTIVE (clean slab): metal atoms within 3.6 A of the point 2.3 A straight above the site atom ("crowding" of the pocket)
    q = P[i] + np.array([0.0, 0.0, 2.3]); dq = P - q
    dq[:, 0] -= Lx * np.round(dq[:, 0] / Lx); dq[:, 1] -= Ly * np.round(dq[:, 1] / Ly)
    out["crowd"] = int((np.sqrt((dq ** 2).sum(1)) < 3.6).sum())
    # EXPLANATORY (relaxed H2O*): metal atoms within 3.0 A of the O, and O height above the top metal atom
    fa = final_atoms(f"qe_inputs/ads_{cell}_{site}_H2O.out")
    o = [a for a in fa if a[0] == "O"]
    if o:
        met = np.array([[a[1], a[2], a[3]] for a in fa if a[0] in ("Mg", "Zn", "Ca")])
        do = met - np.array(o[0][1:]); do[:, 0] -= Lx * np.round(do[:, 0] / Lx); do[:, 1] -= Ly * np.round(do[:, 1] / Ly)
        out["contacts"] = int((np.sqrt((do ** 2).sum(1)) < 3.0).sum()); out["zO"] = float(o[0][3] - met[:, 2].max())
    else:
        out["contacts"] = -1; out["zO"] = float("nan")
    return out


rows = []
for r in csv.DictReader(open("ads/ads_table.csv")):
    if r["mol"] != "H2O" or r["cut"] != "std":
        continue
    try:
        em = float(r["Em"])
    except ValueError:
        continue
    if not np.isfinite(em):
        continue
    sites = {s["site"]: s for s in csv.DictReader(open(f"ads/sites_{r['cell']}.csv"))}
    s = sites[r["site"]]
    d = descriptors(r["cell"], int(s["atom"]), r["site"])
    d.update(cell=r["cell"], site=r["site"], el=s["el"], Em=em, dZ=float(r["dZ"]) if r["dZ"] not in ("", "nan") else np.nan,
             crystal=r["cell"].startswith("crystal"))
    rows.append(d)
cry = [r for r in rows if r["crystal"]]
am = [r for r in rows if not r["crystal"] and np.isfinite(r["dZ"])]
print(f"{len(am)} amorphous sites, {len(cry)} crystal site(s)\n")
print(f"{'cell':14s} {'site':4s} el  CN3.2 CN3.4 CN3.6 depth/A nZnCa crowd cont3.0 zO/A  dZ      E_matched/eV")
for r in sorted(am + cry, key=lambda r: (r["el"], r["cn34"])):
    print(f"{r['cell']:14s} {r['site']:4s} {r['el']}  {r['cn32']:4d} {r['cn34']:5d} {r['cn36']:5d}  {r['depth']:6.2f} {r['nZnCa']:5d} {r['crowd']:5d} {r['contacts']:6d} {r['zO']:5.2f}  {r['dZ']:+.3f}  {r['Em']:+.3f}")


def rank(a):
    return np.argsort(np.argsort(a)).astype(float)


def spearman_p(x, y, n=20000):
    x, y = rank(np.asarray(x)), rank(np.asarray(y))
    rho = np.corrcoef(x, y)[0, 1]
    c = sum(abs(np.corrcoef(x, rng.permutation(y))[0, 1]) >= abs(rho) - 1e-12 for _ in range(n))
    return rho, (c + 1) / (n + 1)


print("\nSpearman rho of E_matched with CN(3.4 A) [permutation p]; other cutoffs in brackets")
for lab, sel in (("all amorphous", am),) + tuple((e, [r for r in am if r["el"] == e]) for e in ELS):
    if len(sel) < 5:
        continue
    y = [r["Em"] for r in sel]
    o = {k: spearman_p([r[k] for r in sel], y) for k in ("cn34", "cn32", "cn36", "depth", "nZnCa", "crowd", "contacts")}
    print(f"  {lab:14s} n={len(sel):2d}  CN {o['cn34'][0]:+.2f} [p={o['cn34'][1]:.3f}]  (3.2: {o['cn32'][0]:+.2f}, 3.6: {o['cn36'][0]:+.2f})"
          f"   depth {o['depth'][0]:+.2f} [p={o['depth'][1]:.3f}]   nZnCa {o['nZnCa'][0]:+.2f} [p={o['nZnCa'][1]:.3f}]"
          f"\n{'':20s} crowd(clean, predictive) {o['crowd'][0]:+.2f} [p={o['crowd'][1]:.3f}]   O-metal contacts<3.0 A (relaxed, explanatory) {o['contacts'][0]:+.2f} [p={o['contacts'][1]:.3f}]")
print("  (E_ads is negative for binding, so a POSITIVE rho with CN means more neighbours -> weaker binding)")

def design(sel, terms):
    cols = [np.array([r["el"] == e for r in sel], float) for e in ELS if any(r["el"] == e for r in am)]
    for t in terms:
        if t.startswith("x:"):                       # element-specific slope
            k = t[2:]
            for e in ELS:
                v = np.array([r[k] if r["el"] == e else 0.0 for r in sel], float)
                if np.any(v):
                    cols.append(v - np.where(v != 0, np.mean([r[k] for r in sel if r["el"] == e]), 0.0))
        else:
            cols.append(np.array([r[t] for r in sel], float) - np.mean([r[t] for r in am]))
    return np.column_stack(cols)


def loo(terms):
    y = np.array([r["Em"] for r in am]); X = design(am, terms); err = []
    for i in range(len(y)):
        m = np.arange(len(y)) != i
        b = np.linalg.lstsq(X[m], y[m], rcond=None)[0]
        err.append(y[i] - X[i] @ b)
    return float(np.sqrt(np.mean(np.square(err))))


print("\nLeave-one-site-out RMSE of E_matched (eV), lower = better")
for lab, terms in (("element only", []), ("element + CN3.4", ["cn34"]), ("element + depth", ["depth"]), ("element + nZnCa", ["nZnCa"]),
                   ("element + dZ", ["dZ"]), ("element + crowd (clean slab)", ["crowd"]), ("element + O contacts (relaxed; explanatory)", ["contacts"]), ("element + CN3.4 + dZ", ["cn34", "dZ"]), ("element x CN3.4 (own slopes)", ["x:cn34"])):
    print(f"  {lab:30s} {loo(terms):.3f}")

print("\nTransfer to the Mg(0001) atop site (DFT E_matched shown first)")
for r in cry:
    print(f"  {r['cell']} {r['site']}: DFT {r['Em']:+.3f} eV, CN3.4 = {r['cn34']}, depth {r['depth']:.2f} A")
    y = np.array([q["Em"] for q in am])
    for lab, terms in (("element only", []), ("element + CN3.4", ["cn34"]), ("element + depth", ["depth"]), ("element + crowd", ["crowd"])):
        X = design(am, terms); b = np.linalg.lstsq(X, y, rcond=None)[0]
        xr = [float(r["el"] == e) for e in ELS if any(q["el"] == e for q in am)]
        xr += [float(r[t]) - np.mean([q[t] for q in am]) for t in terms]
        print(f"    predicted by {lab:18s} {float(np.array(xr) @ b):+.3f} eV")
    mg = [q for q in am if q["el"] == "Mg"]
    if len(mg) > 4:
        for k in ("cn34", "depth"):
            c = np.polyfit([q[k] for q in mg], [q["Em"] for q in mg], 1)
            print(f"    Mg-only fit on {k}: {np.polyval(c, r[k]):+.3f} eV (slope {c[0]:+.3f} eV per unit)")
print("\nMg amorphous CN range:", sorted({q["cn34"] for q in am if q["el"] == "Mg"}), " crystal atop CN:", [r["cn34"] for r in cry])
