"""
Barrier estimate for H2O* -> OH* + H* by a RESTRAINED relaxed scan (pw.x only; no neb.x needed).

  python3 13_barrier.py build     # inputs bar_<cell>_<site>_L<k>.in for one site per element
  python3 13_barrier.py analyze   # energy profile, barrier, checks

Method
  reactant = dref_<cell>_<site> (relaxed molecular water), product = ddis_<cell>_<site> (relaxed OH* + H*), same free atoms.
  Images: positions interpolated linearly between reactant and product (lambda = 0.15 ... 0.90). In each image the MOVING H is
  frozen at its interpolated position (if_pos 0 0 0) and all other free atoms (O, the other H, the nearby slab) relax with BFGS.
  (QE's CONSTRAINTS card is not usable here: "constraints only with fixed-cell dynamics".)
  E_profile(lambda) - E(reactant) maximum = barrier ESTIMATE. A relaxed one-coordinate scan is an UPPER bound of the true
  barrier (freezing the H position is stricter than a distance constraint, and the best path can need other coordinates); good for ranking Ca / Mg / Zn, not for absolute rate constants.
  The profile is only trustworthy if it is smooth, every image converged, and its maximum is interior.
Env: ELEMS ("Ca Mg Zn"), BARRIER_SITES ("cell:site cell:site ..." overrides the automatic choice), LAMBDAS.
Site choice: per element, among dissociation sites with both ends converged, the one whose dE_diss is closest to the element median.
"""
import csv, glob, importlib, os, re, sys
import numpy as np

A = importlib.import_module("11_adsorption")
RY = A.RY
ELEMS = os.environ.get("ELEMS", "Ca Mg Zn").split()
LAMBDAS = [float(x) for x in os.environ.get("LAMBDAS", "0.15 0.30 0.45 0.60 0.75 0.90").split()]


def finished(name, dis=False):
    fn = f"qe_inputs/{name}.out"
    if not os.path.exists(fn):
        return False
    t = open(fn, errors="ignore").read()
    return "JOB DONE" in t and "convergence NOT achieved" not in t and A.last_energy(fn) is not None


def site_element(cell, site):
    fn = f"ads/sites_{cell}.csv"
    for r in csv.DictReader(open(fn)):
        if r["site"] == site:
            return r["el"]
    return None


def candidates():
    out = []
    for fn in sorted(glob.glob("qe_inputs/ddis_*.in")):
        base = os.path.basename(fn)[5:-3]
        m = re.match(r"(.+)_(s\d+)$", base)
        if not m:
            continue
        cell, site = m.groups()
        if finished(f"dref_{base}") and finished(f"ddis_{base}"):
            dE = (A.last_energy(f"qe_inputs/ddis_{base}.out") - A.last_energy(f"qe_inputs/dref_{base}.out")) * RY
            out.append((site_element(cell, site), cell, site, dE))
    return out


def choose():
    env = os.environ.get("BARRIER_SITES")
    cand = candidates()
    print("valid dissociation sites (both ends converged):")
    for el, c, s, d in cand:
        print(f"  {el}  {c}:{s}  dE_diss = {d:+.3f} eV")
    chosen = {}
    if env:
        for item in env.split():
            c, s = item.split(":")
            chosen[site_element(c, s)] = (c, s)
        return chosen
    for el in ELEMS:
        cs = [x for x in cand if x[0] == el]
        if not cs:
            print(f"no valid dissociation site for {el} yet - skipped"); continue
        med = np.median([x[3] for x in cs])
        _, c, s, _ = min(cs, key=lambda x: abs(x[3] - med))
        chosen[el] = (c, s)
    return chosen


def mic(d, L):
    return d - L * np.round(d / L)


def build():
    chosen = choose()
    names = {}
    for el, (cell, site) in chosen.items():
        base = f"{cell}_{site}"
        tmpl = A.load(f"qe_inputs/{cell}_wf.in"); n = len(tmpl["atoms"])
        Lx, Ly = tmpl["cell"][0][0], tmpl["cell"][1][1]
        R = A.final_positions(f"qe_inputs/dref_{base}.out"); Pd = A.final_positions(f"qe_inputs/ddis_{base}.out")
        assert len(R) == n + 3 and len(Pd) == n + 3, "unexpected atom count"
        free = [a[4] for a in A.parse_flags(f"qe_inputs/dref_{base}.in")]
        # product order: slab, O, H(stay), H(moved). Reactant H that is closer to product H(stay) is the one that stays.
        def d3(a, b):
            return np.sqrt(mic(a[1] - b[1], Lx) ** 2 + mic(a[2] - b[2], Ly) ** 2 + (a[3] - b[3]) ** 2)
        h1, h2 = R[n + 1], R[n + 2]
        if d3(h1, Pd[n + 1]) <= d3(h2, Pd[n + 1]):
            Rr = R[:n + 1] + [h1, h2]
        else:
            Rr = R[:n + 1] + [h2, h1]
        names[el] = []
        for k, lam in enumerate(LAMBDAS, 1):
            name = f"bar_{base}_L{k}"
            if finished(name):
                names[el].append(name); continue
            atoms = []
            for a, b in zip(Rr, Pd):
                dx, dy, dz = mic(b[1] - a[1], Lx), mic(b[2] - a[2], Ly), b[3] - a[3]
                atoms.append((a[0], a[1] + lam * dx, a[2] + lam * dy, a[3] + lam * dz))
            fimg = list(free); fimg[n + 2] = False                      # the moving H is frozen at its interpolated position
            A.write_input(f"qe_inputs/{name}.in", tmpl, atoms, fimg, name, "relax", extra_bands=A.EXTRA_BANDS, new_species=("H", "O"))
            t = open(f"qe_inputs/{name}.in").read()
            t = re.sub(r"mixing_beta *= *[0-9.]+", "mixing_beta = 0.1", t)
            t = re.sub(r"\n *electron_maxstep *= *\d+", "", t)
            t = t.replace("&electrons\n", "&electrons\n  electron_maxstep = 300\n", 1)
            open(f"qe_inputs/{name}.in", "w").write(t)
            names[el].append(name)
        print(f"{el}: {cell}:{site}  ->  {len(LAMBDAS)} constrained images, {sum(free)} free atoms")
    for el, lst in names.items():
        open(f"ads/list_bar_{el}.txt", "w").write("\n".join(lst) + "\n")
    print("lists:", {el: len(v) for el, v in names.items()})


REFINE = {"Mg": [0.22, 0.26, 0.34, 0.38, 0.17, 0.19], "Zn": [0.22, 0.26, 0.34, 0.38, 0.36],
          "Ca": [0.38, 0.42, 0.50, 0.54, 0.33, 0.35]}   # round 1 around each maximum; round 2 (last entries) brackets the edge maxima


def refine():
    """Extra images near the profile maximum, named bar_<cell>_<site>_P<100*lambda>. Each starts from the INTERPOLATION OF THE
    TWO NEIGHBOURING RELAXED IMAGES (not from the raw reactant/product line), which keeps the path continuous. Moving H frozen again.
    Same site as the main scan (choose()); override lambdas with REFINE_LAMBDAS="0.2 0.3" (applies to all elements)."""
    chosen = choose()
    names = {}
    for el, (cell, site) in chosen.items():
        base = f"{cell}_{site}"
        tmpl = A.load(f"qe_inputs/{cell}_wf.in"); n = len(tmpl["atoms"])
        Lx, Ly = tmpl["cell"][0][0], tmpl["cell"][1][1]
        R = A.final_positions(f"qe_inputs/dref_{base}.out"); Pd = A.final_positions(f"qe_inputs/ddis_{base}.out")
        def d3(a, b):
            return np.sqrt(mic(a[1] - b[1], Lx) ** 2 + mic(a[2] - b[2], Ly) ** 2 + (a[3] - b[3]) ** 2)
        h1, h2 = R[n + 1], R[n + 2]
        Rr = R[:n + 1] + ([h1, h2] if d3(h1, Pd[n + 1]) <= d3(h2, Pd[n + 1]) else [h2, h1])
        known = {0.0: Rr, 1.0: Pd}
        for k, lam in enumerate(LAMBDAS, 1):
            out = f"qe_inputs/bar_{base}_L{k}.out"
            if finished(f"bar_{base}_L{k}") and "bfgs converged" in open(out, errors="ignore").read():
                known[lam] = A.final_positions(out)
        for fp in glob.glob(f"qe_inputs/bar_{base}_P*.out"):               # finished refined images are neighbours too
            pm = re.search(r"_P(\d+)\.out$", fp)
            if pm and "bfgs converged" in open(fp, errors="ignore").read() and "JOB DONE" in open(fp, errors="ignore").read():
                known[int(pm.group(1)) / 100.0] = A.final_positions(fp)
        free = [a[4] for a in A.parse_flags(f"qe_inputs/dref_{base}.in")]
        fimg = list(free); fimg[n + 2] = False
        lams = [float(x) for x in os.environ["REFINE_LAMBDAS"].split()] if os.environ.get("REFINE_LAMBDAS") else REFINE.get(el, [])
        names[el] = []
        for lam in lams:
            name = f"bar_{base}_P{int(round(lam * 100))}"
            if finished(name):
                names[el].append(name); continue
            lo = max(l for l in known if l < lam); hi = min(l for l in known if l > lam)
            t = (lam - lo) / (hi - lo)
            atoms = []
            for a, b in zip(known[lo], known[hi]):
                dx, dy, dz = mic(b[1] - a[1], Lx), mic(b[2] - a[2], Ly), b[3] - a[3]
                atoms.append((a[0], a[1] + t * dx, a[2] + t * dy, a[3] + t * dz))
            A.write_input(f"qe_inputs/{name}.in", tmpl, atoms, fimg, name, "relax", extra_bands=A.EXTRA_BANDS, new_species=("H", "O"))
            tx = open(f"qe_inputs/{name}.in").read()
            tx = re.sub(r"mixing_beta *= *[0-9.]+", "mixing_beta = 0.1", tx)
            tx = re.sub(r"\n *electron_maxstep *= *\d+", "", tx)
            tx = tx.replace("&electrons\n", "&electrons\n  electron_maxstep = 300\n", 1)
            tx = re.sub(r"nstep *= *\d+", "nstep = 200", tx)
            open(f"qe_inputs/{name}.in", "w").write(tx)
            names[el].append(name)
            print(f"{el}: {name}  (between lambda {lo:.2f} and {hi:.2f})")
    for el, lst in names.items():
        open(f"ads/list_bar_ref_{el}.txt", "w").write("\n".join(lst) + ("\n" if lst else ""))
    print("lists:", {el: len(v) for el, v in names.items()})


def cont():
    """Restart scan images that finished their 80 BFGS steps without converging: restart from the last geometry with nstep=200.
    The unconverged output is kept as <name>.out.part1 (appended parts: part2, ...). Then run them with submit_ads.sh (list ads/list_bar_cont.txt)."""
    todo = []
    for fn in sorted(glob.glob("qe_inputs/bar_*.out")):
        name = os.path.basename(fn)[:-4]
        t = open(fn, errors="ignore").read()
        if "JOB DONE" in t and "bfgs converged" not in t:
            todo.append(name)
    for name in todo:
        fin = A.final_positions(f"qe_inputs/{name}.out")
        old = A.parse_flags(f"qe_inputs/{name}.in")
        assert len(fin) == len(old), name
        txt = open(f"qe_inputs/{name}.in").read()
        head, rest = txt.split("ATOMIC_POSITIONS angstrom", 1)
        _, tail = rest.split("K_POINTS", 1)
        pos = "".join("  %s  %.6f  %.6f  %.6f%s\n" % (a[0], a[1], a[2], a[3], "  1 1 1" if o[4] else "  0 0 0") for a, o in zip(fin, old))
        head = re.sub(r"nstep *= *\d+", "nstep = 200", head)
        open(f"qe_inputs/{name}.in", "w").write(head + "ATOMIC_POSITIONS angstrom\n" + pos + "\nK_POINTS" + tail)
        k = 1
        while os.path.exists(f"qe_inputs/{name}.out.part{k}"):
            k += 1
        os.rename(f"qe_inputs/{name}.out", f"qe_inputs/{name}.out.part{k}")
        print(f"restart {name} from its last geometry (part{k} kept)")
    if todo:                       # never overwrite the list with an empty one (a second call finds nothing: the outputs were moved)
        open("ads/list_bar_cont.txt", "w").write("\n".join(todo) + "\n")
    print(f"{len(todo)} images to continue")


def analyze():
    found = False
    for fn in sorted(glob.glob("qe_inputs/dref_*.in")):
        base = os.path.basename(fn)[5:-3]
        imgs = []
        for im in glob.glob(f"qe_inputs/bar_{base}_[LP]*.in"):
            mm = re.search(r"_([LP])(\d+)\.in$", im)
            lam = LAMBDAS[int(mm.group(2)) - 1] if mm.group(1) == "L" else int(mm.group(2)) / 100.0
            imgs.append((lam, im))
        imgs.sort()
        if not imgs or not finished(f"dref_{base}"):
            continue
        found = True
        m = re.match(r"(.+)_(s\d+)$", base); cell, site = m.groups()
        el = site_element(cell, site)
        e0 = A.last_energy(f"qe_inputs/dref_{base}.out")
        pts = [(0.0, 0.0, "reactant", True)]
        for lam, im in imgs:
            name = os.path.basename(im)[:-3]; out = f"qe_inputs/{name}.out"
            if not os.path.exists(out) or A.last_energy(out) is None:
                pts.append((lam, float("nan"), "pending", False)); continue
            t = open(out, errors="ignore").read()
            ok = ("JOB DONE" in t) and ("convergence NOT achieved" not in t) and ("bfgs converged" in t)
            pts.append((lam, (A.last_energy(out) - e0) * RY, "ok" if ok else "NOT CONVERGED", ok))
        if finished(f"ddis_{base}"):
            pts.append((1.0, (A.last_energy(f"qe_inputs/ddis_{base}.out") - e0) * RY, "product", True))
        print(f"\n{el}  {cell}:{site}   E relative to molecular H2O* (eV)")
        for lam, e, st, ok in pts:
            print(f"  lambda={lam:4.2f}  {e:+7.3f}  {st}")
        good = [(l, e) for l, e, s, ok in pts if ok and np.isfinite(e)]
        con = [(l, e) for l, e in good if 0 < l < 1]
        if con:
            lm, bm = max(con, key=lambda x: x[1])
            interior = 0 < [l for l, _ in good].index(lm) < len(good) - 1 and lm != con[-1][0] if len(con) > 1 else False
            dE = [e for l, e in good if l == 1.0]
            msg = f"  highest sampled image = {max(bm, 0):.3f} eV at lambda={lm:.2f} (barrier estimate; not a strict bound)"
            if dE:
                lowest = min(e for _, e in good)
                msg += f"; relaxed product {dE[0]:+.3f} eV, lowest image {lowest:+.3f} eV"
                if lowest < dE[0] - 0.05:
                    msg += " [an image lies below the relaxed product: the product is a shallow local minimum]"
            if lm == con[-1][0] and len(con) > 1:
                msg += "  [maximum at the last image: profile not past the top - add images closer to lambda=1]"
            if any(st == "pending" for _, _, st, _ in pts):
                msg += "  [incomplete profile]"
            print(msg)
    if not found:
        print("no barrier scans finished yet")


if __name__ == "__main__":
    {"build": build, "analyze": analyze, "cont": cont, "refine": refine}.get(sys.argv[1] if len(sys.argv) > 1 else "", lambda: print(__doc__))()
