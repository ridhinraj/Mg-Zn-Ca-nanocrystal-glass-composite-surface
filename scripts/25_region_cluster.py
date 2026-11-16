"""
Is the OH interphase-vs-glass difference real, or pseudo-replication?  (25_region_cluster.py)

11_adsorption.py analyze reports, at site level:
    OH   interphase minus glass matrix  -0.329 eV  p = 0.012  (19 vs 15 sites)
    H2O  interphase minus glass matrix  -0.013 eV  p = 0.881
A significant region effect for OH would cut against the paper's central null. But that test treats
sites as independent units, and they are not: the six interphase patches are three cut depths at each
of two poles, sharing glass atoms, so the manuscript already calls them "three independent realisations
at two poles, not six replicates". Several sites also come from one patch.

This recomputes the comparison at three levels of aggregation, permuting whole CLUSTERS rather than
sites, and checks the coordination confound (the interphase offers more 3-fold hollows than the glass):

  site   every adsorption site its own unit        - anti-conservative, what analyze prints
  patch  one value per patch (mean of its sites)   - removes within-patch sharing
  pole   iface_N / iface_S / each matrix patch     - the manuscript's own independent realisations

  python3 25_region_cluster.py
"""
import csv, os, re
import numpy as np

rng = np.random.default_rng(0)
F = lambda v: float(v) if str(v).strip() not in ("", "nan", "None") else float("nan")


def region(c):
    return "interphase" if c.startswith("iface_") else "glass matrix" if c.startswith("matrix_") else "crystal"


def pole(c):
    m = re.match(r"(iface_[NS])_D\d+", c)
    return m.group(1) if m else c


def load():
    out = []
    for r in csv.DictReader(open("ads/ads_table.csv")):
        if r["cut"] != "std" or region(r["cell"]) == "crystal":
            continue
        e = F(r["Em"])
        if np.isfinite(e):
            out.append({"cell": r["cell"], "site": r["site"], "mol": r["mol"], "el": r["el"],
                        "Em": e, "cn": F(r.get("cn"))})
    return out


def cluster_perm(vals, labels, n=20000):
    """permute whole clusters between the two regions"""
    keys = sorted(vals)
    a = [vals[k] for k in keys if labels[k] == "interphase"]
    b = [vals[k] for k in keys if labels[k] == "glass matrix"]
    if len(a) < 2 or len(b) < 2:
        return np.mean(a) - np.mean(b) if a and b else float("nan"), float("nan"), len(a), len(b)
    obs = np.mean(a) - np.mean(b)
    pool = np.array([vals[k] for k in keys]); k = len(a)
    c = sum(abs((q := rng.permutation(pool))[:k].mean() - q[k:].mean()) >= abs(obs) - 1e-12 for _ in range(n))
    return obs, (c + 1) / (n + 1), len(a), len(b)


D = load()
if not D:
    raise SystemExit("no amorphous rows in ads/ads_table.csv")

for mol in ("H2O", "OH"):
    v = [r for r in D if r["mol"] == mol]
    if len(v) < 6:
        continue
    print("=" * 78)
    print(f"{mol}   {len(v)} amorphous sites")
    print("=" * 78)
    for lab, keyf in (("site", lambda r: f"{r['cell']}:{r['site']}"),
                      ("patch", lambda r: r["cell"]),
                      ("pole", lambda r: pole(r["cell"]))):
        g = {}
        for r in v:
            g.setdefault(keyf(r), []).append(r["Em"])
        vals = {k: float(np.mean(a)) for k, a in g.items()}
        labels = {k: region(k.split(":")[0] if lab == "site" else k if lab == "patch"
                            else (k if k.startswith("matrix") else k + "_D0")) for k in g}
        d, p, na, nb = cluster_perm(vals, labels)
        print(f"  {lab:6s} units {na:3d} interphase vs {nb:3d} glass   "
              f"difference {d:+.3f} eV   cluster-permutation p = "
              + (f"{p:.4f}" if np.isfinite(p) else "n/a"))
    # coordination confound
    print("  coordination mix (share of that region's sites)")
    for r_ in ("interphase", "glass matrix"):
        s = [x for x in v if region(x["cell"]) == r_]
        cc = {}
        for x in s:
            if np.isfinite(x["cn"]):
                cc[int(x["cn"])] = cc.get(int(x["cn"]), 0) + 1
        print(f"    {r_:13s} " + "  ".join(f"cn={c}: {100*n/len(s):.0f}% ({n})" for c, n in sorted(cc.items())))
    print("  region difference at fixed coordination (sites, indicative only)")
    for c in sorted({int(x["cn"]) for x in v if np.isfinite(x["cn"])}):
        a = [x["Em"] for x in v if np.isfinite(x["cn"]) and int(x["cn"]) == c and region(x["cell"]) == "interphase"]
        b = [x["Em"] for x in v if np.isfinite(x["cn"]) and int(x["cn"]) == c and region(x["cell"]) == "glass matrix"]
        if len(a) >= 2 and len(b) >= 2:
            print(f"    cn={c}: interphase {np.mean(a):+.3f} (n={len(a)})  glass {np.mean(b):+.3f} (n={len(b)})"
                  f"   difference {np.mean(a) - np.mean(b):+.3f}")
    print()

print("""Reading it:
  If the difference shrinks sharply from site to patch to pole and loses significance, the site-level
  result is pseudo-replication and the region null stands for OH as it does for H2O - but the
  manuscript should then say the region comparison is made on independent realisations, not on sites,
  because a referee who reruns the site-level test will get p = 0.012 and ask.
  If it survives aggregation, OH does distinguish the regions even though H2O does not, and that
  belongs in the paper: the interphase would matter for the hydroxide-covered surface but not for
  initiation - a qualification of the title claim, not a refutation of it.""")
