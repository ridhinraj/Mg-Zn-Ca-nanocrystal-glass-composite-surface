"""
OH at every existing H2O site - element-stratified by construction, and paired. (21_oh.py)

OH rather than H2O is the species on a corroding Mg surface at physiological pH, and the element effect is
established for H2O and merely untested for OH. This puts OH on the SAME sites the H2O calibration used, so
the two are paired site by site and the element strata carry over unchanged.

It is strictly ADDITIVE and never re-picks or renumbers anything:
  * sites come from the existing ads/sites_<cell>.csv - 11_adsorption.py build is NOT re-run
  * clean references reuse an existing cref whose free-atom set matches; a new one is minted as
    cref_<cell>_oh<k>, so no existing cref_<cell>_<n>.in is ever overwritten
  * 11_adsorption.py cleanref must NOT be re-run on these cells: it names crefs by position in a sorted
    glob, so adding OH inputs would renumber them and silently invalidate every finished H2O energy

  python3 21_oh.py plan          # what exists, what would be built - writes nothing
  python3 21_oh.py build         # OH inputs + any new crefs + job lists
  python3 21_oh.py report        # paired H2O vs OH, by element

  CELLS="..."   limit to these cells (default: every ads/sites_*.csv)
  MAX_SITES=N   cap the number of new OH sites, keeping the element strata balanced
  NJOBS=N       how many job lists to split across (default 3)
"""
import csv, glob, importlib.util, os, re, sys
import numpy as np

spec = importlib.util.spec_from_file_location("ads11", os.path.join(os.path.dirname(os.path.abspath(__file__)), "11_adsorption.py"))
A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)

ELS = ("Ca", "Mg", "Zn")
NJOBS = int(os.environ.get("NJOBS", "3"))
MAX_SITES = int(os.environ.get("MAX_SITES", "0"))
rng = np.random.default_rng(0)


def cells():
    c = os.environ.get("CELLS", "").split()
    return c or sorted(os.path.basename(f)[6:-4] for f in glob.glob("ads/sites_*.csv"))


def sites_of(cell):
    fn = f"ads/sites_{cell}.csv"
    return list(csv.DictReader(open(fn))) if os.path.exists(fn) else []


def done(name):
    fn = f"qe_inputs/{name}.out"
    return os.path.exists(fn) and "JOB DONE" in open(fn, errors="ignore").read()


def slab_key(fn, n):
    return tuple(a[4] for a in A.parse_flags(fn)[:n])


def survey():
    """per cell: template, n, existing cref free-sets, and each site's H2O/OH state"""
    out = []
    for cell in cells():
        ref = f"qe_inputs/{cell}_wf.in"
        if not os.path.exists(ref):
            print(f"skip {cell}: no {ref}"); continue
        tmpl = A.load(ref); n = len(tmpl["atoms"])
        crefs = {}
        for cf in sorted(glob.glob(f"qe_inputs/cref_{cell}_*.in")):
            cn = os.path.basename(cf)[:-3]
            try:
                crefs[slab_key(cf, n)] = cn
            except Exception:
                pass
        rows = []
        for s in sites_of(cell):
            s = dict(s); s["atom"] = int(s["atom"]) - 1
            h2o, oh = f"ads_{cell}_{s['site']}_H2O", f"ads_{cell}_{s['site']}_OH"
            rows.append((s, h2o, oh, done(h2o), os.path.exists(f"qe_inputs/{oh}.in"), done(oh)))
        out.append((cell, tmpl, n, crefs, rows))
    return out


def plan(verbose=True):
    sv = survey()
    todo = []
    nh = no = 0
    for cell, tmpl, n, crefs, rows in sv:
        for s, h2o, oh, h2o_done, oh_in, oh_done in rows:
            nh += h2o_done; no += oh_done
            if h2o_done and not oh_done:
                todo.append((cell, tmpl, n, crefs, s, oh, oh_in))
    if MAX_SITES and len(todo) > MAX_SITES:                 # keep the strata balanced
        byel = {}
        for t in todo:
            byel.setdefault(t[4]["el"], []).append(t)
        for k in byel:
            byel[k] = [byel[k][j] for j in rng.permutation(len(byel[k]))]
        keep, i = [], 0
        while len(keep) < MAX_SITES and any(byel.values()):
            for e in ELS:
                if byel.get(e) and len(keep) < MAX_SITES:
                    keep.append(byel[e].pop(0))
        todo = keep
    if verbose:
        print(f"{len(sv)} cells | H2O finished {nh} | OH finished {no} | OH to build {len(todo)}\n")
        for cell, tmpl, n, crefs, rows in sv:
            done_h = sum(1 for r in rows if r[3]); done_o = sum(1 for r in rows if r[5])
            pend = [r[0]["site"] for r in rows if r[3] and not r[5]]
            print(f"  {cell:16s} {len(rows):2d} sites | H2O {done_h:2d} | OH {done_o:2d} | "
                  f"{len(crefs):2d} crefs | to build: {','.join(pend) if pend else '-'}")
        c = {}
        for t in todo:
            c[t[4]["el"]] = c.get(t[4]["el"], 0) + 1
        print("\nnew OH sites by element: " + (", ".join(f"{e} {c.get(e,0)}" for e in ELS) if todo else "none"))
        skipped = [(cell, r[0]["site"]) for cell, _, _, _, rows in sv for r in rows if not r[3]]
        if skipped:
            print(f"\n{len(skipped)} sites have no finished H2O run and are left alone: "
                  + ", ".join(f"{c}:{s}" for c, s in skipped[:12]) + (" ..." if len(skipped) > 12 else ""))
    return todo


def build():
    todo = plan(verbose=True)
    if not todo:
        print("\nnothing to build"); return
    newc, names, crefnames = {}, [], []
    for cell, tmpl, n, crefs, s, oh, oh_in in todo:
        atoms, free = A.place(tmpl, s, "OH")
        A.write_input(f"qe_inputs/{oh}.in", tmpl, atoms, free, oh, "relax",
                      extra_bands=A.EXTRA_BANDS, new_species=("H", "O"))
        names.append(oh)
        key = tuple(free[:n])
        pool = dict(crefs); pool.update(newc.get(cell, {}))
        if key in pool:
            cn = pool[key]
        else:
            k = len(newc.setdefault(cell, {})) + 1
            while os.path.exists(f"qe_inputs/cref_{cell}_oh{k}.in"):
                k += 1
            cn = f"cref_{cell}_oh{k}"
            A.write_input(f"qe_inputs/{cn}.in", tmpl, tmpl["atoms"], list(key), cn, "relax")
            newc[cell][key] = cn
            crefnames.append(cn)
        with open("ads/cref_map.csv", "a", newline="") as f:
            csv.writer(f).writerow([oh, cn])
    allj = names + crefnames
    for k in range(NJOBS):
        part = allj[k::NJOBS]
        open(f"ads/list_oh_{k+1}.txt", "w").write("\n".join(part) + ("\n" if part else ""))
    print(f"\n{len(names)} OH inputs, {len(crefnames)} new clean references "
          f"(reused {len(names) - len(crefnames)} existing)")
    for k in range(NJOBS):
        print(f"  ads/list_oh_{k+1}.txt: {len(allj[k::NJOBS])} jobs")
    print("\nappended to ads/cref_map.csv; no existing cref was touched.\nnext:  bash run_oh.sh")


def report():
    cmap = {}
    if os.path.exists("ads/cref_map.csv"):
        for r in csv.DictReader(open("ads/cref_map.csv")):
            cmap[r["ads"]] = r["cref"]                      # later rows win
    def E(nm):
        fn = f"qe_inputs/{nm}.out"
        return A.last_energy(fn) * A.RY if os.path.exists(fn) and done(nm) else None
    gw, gh = E("gas_H2O"), E("gas_H2")
    if gw is None or gh is None:
        print("gas_H2O / gas_H2 not finished - cannot reference"); return
    ref = {"H2O": gw, "OH": gw - 0.5 * gh}
    rows = []
    for cell, tmpl, n, crefs, srows in survey():
        for s, h2o, oh, h2o_done, oh_in, oh_done in srows:
            rec = {"cell": cell, "site": s["site"], "el": s["el"]}
            for mol, nm in (("H2O", h2o), ("OH", oh)):
                e, c = E(nm), cmap.get(nm)
                ec = E(c) if c else None
                rec[mol] = (e - ec - ref[mol]) if (e is not None and ec is not None) else None
            rows.append(rec)
    # The model is applied to amorphous surfaces only, so crystal_ideal is reported separately and
    # kept OUT of the element means. Its hollow sites bind OH at -1.62 eV against -0.30 eV atop, so
    # including them moves the Mg mean by ~0.06 eV and inflates its scatter.
    ok = [r for r in rows if r["OH"] is not None and not r["cell"].startswith("crystal")]
    cry = [r for r in rows if r["OH"] is not None and r["cell"].startswith("crystal")]
    print(f"{len(ok)} amorphous OH sites finished "
          f"({sum(1 for r in ok if r['H2O'] is not None)} with H2O at the same site)"
          + (f", plus {len(cry)} crystal sites shown separately" if cry else "") + "\n")
    print(f"{'cell':16s} {'site':5s} {'el':3s} {'E_ads(H2O)':>11s} {'E_ads(OH)':>10s} {'OH-H2O':>8s}")
    for r in sorted(ok + cry, key=lambda r: (r["cell"].startswith("crystal"), r["el"], r["cell"], r["site"])):
        d = (r["OH"] - r["H2O"]) if r["H2O"] is not None else None
        print(f"{r['cell']:16s} {r['site']:5s} {r['el']:3s} "
              + (f"{r['H2O']:+11.3f}" if r["H2O"] is not None else f"{'-':>11s}")
              + f" {r['OH']:+10.3f} " + (f"{d:+8.3f}" if d is not None else f"{'-':>8s}"))

    def perm(d):
        d = np.asarray([x for x in d if x is not None], float)
        if len(d) < 3:
            return np.nan, np.nan, np.nan, len(d)
        m = d.mean()
        c = sum(abs((d * rng.choice([-1, 1], len(d))).mean()) >= abs(m) - 1e-12 for _ in range(20000))
        return m, d.std(ddof=1) / np.sqrt(len(d)), (c + 1) / 20001, len(d)

    def perm2(x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        obs = abs(y.mean() - x.mean()); pool = np.concatenate([x, y]); k = len(x)
        c = sum(abs((q := rng.permutation(pool))[k:].mean() - q[:k].mean()) >= obs - 1e-12 for _ in range(20000))
        return (c + 1) / 20001

    print("\nE_ads(OH) by element (site-matched)")
    byel = {}
    for r in ok:
        byel.setdefault(r["el"], []).append(r["OH"])
    for e in ELS:
        v = byel.get(e, [])
        if v:
            sd = np.std(v, ddof=1) if len(v) > 1 else float("nan")
            print(f"  {e:3s} {np.mean(v):+.3f} +- {sd:.3f} eV   n={len(v)}")
    pairs = [(a, b) for a, b in (("Zn", "Ca"), ("Zn", "Mg"), ("Mg", "Ca"))
             if len(byel.get(a, [])) >= 3 and len(byel.get(b, [])) >= 3]
    for a, b in pairs:
        print(f"  {a} vs {b}: difference {np.mean(byel[b]) - np.mean(byel[a]):+.3f} eV, "
              f"permutation p = {perm2(byel[a], byel[b]):.4f}")
    if not pairs:
        print("  (no element pair has 3 sites on both sides yet)")

    print("\nOH minus H2O at the same site (paired)")
    m, se, p, k = perm([(r["OH"] - r["H2O"]) for r in ok if r["H2O"] is not None])
    print(f"  all sites      {m:+.3f} +- {se:.3f} eV   sign-permutation p = {p:.4f}   (n={k})")
    for e in ELS:
        d = [(r["OH"] - r["H2O"]) for r in ok if r["el"] == e and r["H2O"] is not None]
        m, se, p, k = perm(d)
        if k >= 3:
            print(f"  {e:3s}            {m:+.3f} +- {se:.3f} eV   p = {p:.4f}   (n={k})")
    print("""
What to look for: whether the element ordering Ca < Mg < Zn that holds for H2O also holds for OH, and
whether OH-minus-H2O is a constant offset. A constant offset means the H2O element effect transfers to the
species actually present at physiological pH, which is what the Limitations section currently cannot claim.""")


if __name__ == "__main__":
    os.makedirs("ads", exist_ok=True)
    {"plan": lambda: plan(True), "build": build, "report": report}.get(
        sys.argv[1] if len(sys.argv) > 1 else "", lambda: print(__doc__))()
