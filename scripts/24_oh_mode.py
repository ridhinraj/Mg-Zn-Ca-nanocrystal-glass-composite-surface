"""
Why is OH binding so scattered, and does binding GEOMETRY explain it better than the element?

The OH element means have standard deviations near 0.4 eV against 0.03-0.17 eV for H2O, and the values
are bimodal. The crystal calibration says why it might be: on Mg(0001) OH binds -0.30 eV atop and
-1.62 eV in a hollow. If OH relaxes into whatever hollow the local geometry offers, that 1.3 eV swamps
any element effect - and the right conclusion is not "no element effect for OH" but "for OH the binding
site, not the element under it, sets the energy".

Everything needed is already in ads/ads_table.csv: `cn` is the number of metal atoms coordinating the
adsorbate O after relaxation, `zO` its height. Refresh that file first so the new OH runs are in it:
    sbatch submit_ads_analyze.sh          # runs 11_adsorption.py analyze (read-only; does NOT re-pick
                                          # sites or renumber crefs - that is cleanref, which must not run)

  python3 24_oh_mode.py
"""
import csv, os
import numpy as np

ELS = ("Ca", "Mg", "Zn")
rng = np.random.default_rng(0)
F = lambda v: float(v) if str(v).strip() not in ("", "nan", "None") else float("nan")


def load(amorphous_only=True):
    fn = "ads/ads_table.csv"
    if not os.path.exists(fn):
        raise SystemExit(f"{fn} not found")
    out = []
    for r in csv.DictReader(open(fn)):
        if r["cut"] != "std":
            continue
        if amorphous_only and r["cell"].startswith("crystal"):
            continue
        out.append({"cell": r["cell"], "site": r["site"], "mol": r["mol"], "el": r["el"],
                    "Em": F(r["Em"]), "cn": F(r.get("cn")), "zO": F(r.get("zO"))})
    return [r for r in out if np.isfinite(r["Em"])]


def perm2(x, y, n=20000):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 2 or len(y) < 2:
        return float("nan")
    obs = abs(y.mean() - x.mean()); pool = np.concatenate([x, y]); k = len(x)
    c = sum(abs((q := rng.permutation(pool))[k:].mean() - q[:k].mean()) >= obs - 1e-12 for _ in range(n))
    return (c + 1) / (n + 1)


D = load()
cry = load(amorphous_only=False)
cry = [r for r in cry if r["cell"].startswith("crystal")]
print(f"{len(D)} amorphous rows ({sum(1 for r in D if r['mol']=='OH')} OH, "
      f"{sum(1 for r in D if r['mol']=='H2O')} H2O), {len(cry)} crystal rows\n")
if cry:
    print("crystal_ideal, where geometry is known exactly:")
    for r in sorted(cry, key=lambda r: (r["mol"], r["site"])):
        print(f"   {r['mol']:3s} {r['site']:5s} cn={r['cn']:.0f}  {r['Em']:+.3f} eV")
    print()

for mol in ("H2O", "OH"):
    v = [r for r in D if r["mol"] == mol]
    if len(v) < 4:
        continue
    print(f"=== {mol}  ({len(v)} amorphous sites) " + "=" * 40)
    print("  by element")
    byel = {}
    for r in v:
        byel.setdefault(r["el"], []).append(r["Em"])
    for e in ELS:
        if byel.get(e):
            a = byel[e]
            print(f"    {e:3s} {np.mean(a):+.3f} +- {np.std(a, ddof=1) if len(a)>1 else float('nan'):.3f}  n={len(a)}")
    print("  by coordination of the adsorbate O")
    bycn = {}
    for r in v:
        if np.isfinite(r["cn"]):
            bycn.setdefault(int(r["cn"]), []).append(r["Em"])
    for c in sorted(bycn):
        a = bycn[c]
        print(f"    cn={c}  {np.mean(a):+.3f} +- {np.std(a, ddof=1) if len(a)>1 else float('nan'):.3f}  n={len(a)}"
              + "   " + " ".join(f"{x:+.2f}" for x in sorted(a)))
    # which explains more: element or coordination?
    def rss(groups):
        return sum(float(np.sum((np.array(a) - np.mean(a)) ** 2)) for a in groups.values())
    y = np.array([r["Em"] for r in v]); tot = float(np.sum((y - y.mean()) ** 2))
    if tot > 0:
        print(f"  variance explained:  element {1 - rss(byel)/tot:5.2f}"
              + (f"   coordination {1 - rss(bycn)/tot:5.2f}" if bycn else ""))
    # element effect WITHIN a fixed coordination
    print("  element effect at fixed coordination")
    done = False
    for c in sorted(bycn):
        sub = {}
        for r in v:
            if np.isfinite(r["cn"]) and int(r["cn"]) == c:
                sub.setdefault(r["el"], []).append(r["Em"])
        sub = {k: a for k, a in sub.items() if len(a) >= 2}
        if len(sub) >= 2:
            done = True
            s = "  ".join(f"{k} {np.mean(a):+.3f} (n={len(a)})" for k, a in sorted(sub.items()))
            ks = sorted(sub)
            print(f"    cn={c}: {s}   {ks[0]} vs {ks[-1]} p = {perm2(sub[ks[0]], sub[ks[-1]]):.3f}")
    if not done:
        print("    (no coordination class has two elements with >=2 sites each)")
    print()

print("""Reading it:
  If OH's variance explained by coordination far exceeds that by element, and the element effect
  inside a fixed coordination class is small, then OH binding is set by where the O ends up, not by
  what sits under it. The honest claim then becomes: the exposed element governs the FIRST step
  (molecular H2O adsorption), while for OH the binding geometry dominates - which bounds the model
  to initiation rather than to the hydroxide-covered steady state, and is a result, not a gap.""")
