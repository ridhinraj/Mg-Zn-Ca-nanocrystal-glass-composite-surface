"""
Adsorption calibration (H2O and OH) on selected sites, with Quantum ESPRESSO.

  python3 11_adsorption.py build     # sites + all inputs (cells: CELLS env, default "crystal_ideal matrix_xp")
  python3 11_adsorption.py analyze   # adsorption energies and the descriptor table
  python3 11_adsorption.py dissoc    # H2O* -> OH* + H* inputs for finished sites (DISSOC_SITES)
  python3 11_adsorption.py cleanref  # inputs for site-matched clean-slab references (same free atoms, no adsorbate)

What it builds (all with the SAME cutoffs/smearing/Gamma point as the production runs plus the slab
dipole correction, taken from qe_inputs/<cell>_wf.in):
  * calculation='relax' inputs  qe_inputs/ads_<cell>_<site>_<H2O|OH>.in
        only the adsorbate and the slab atoms within 4.8 A of its O are free; everything else is fixed
  * gas references             qe_inputs/gas_H2O.in, gas_H2.in (14 A box, Gamma)
  * crystal_ideal: sites atop / hcp hollow / fcc hollow; H2O atop, OH on all three
  * amorphous cells: 6 sites picked from the descriptor range (rank of the p-window centre e_p_win among
    protruding surface atoms, >= 4 A apart), H2O and OH atop each
  * a cutoff test on crystal_ideal: H2O atop, clean slab and gas at 60/480 Ry ("_hi")
Energies (eV):  E_ads(H2O) = E(slab+H2O) - E(slab) - E(H2O)
                E_ads(OH)  = E(slab+OH)  - E(slab) - [E(H2O) - 0.5 E(H2)]      (negative = binds)
E(slab) is the dipole-corrected clean SCF in qe_inputs/<cell>_wf.out (MEAM geometry, not DFT-relaxed).
"""
import csv, glob, os, re, sys
import numpy as np

RY = 13.605693123
ZVAL = {"Mg": 2.0, "Zn": 20.0, "Ca": 10.0}   # nominal valence of the SSSP pseudopotentials; dZ = Lowdin q - ZVAL
PSEUDO_FILES = {"H": os.environ.get("PSEUDO_H", "H.pbe-rrkjus_psl.1.0.0.UPF"),
                "O": os.environ.get("PSEUDO_O", "O.pbe-n-kjpaw_psl.0.1.UPF")}
MASS = {"H": 1.008, "O": 15.999}
CELLS = os.environ.get("CELLS", "crystal_ideal matrix_xp").split()
N_SITES = int(os.environ.get("N_SITES", "6"))
FREE_RADIUS = 4.8          # A around the adsorbate O
MIN_SEP = 4.0              # A between selected sites
DEPTH_MAX = float(os.environ.get("DEPTH_MAX", "1.5"))   # only atoms within this depth of the top atom are sites
HI = (60.0, 480.0)         # cutoff test (SSSP recommendation for H)
EXTRA_BANDS = 6
SITE_MODE = os.environ.get("SITE_MODE", "descriptor")      # "element": round-robin Ca, Zn, Mg (random within an element)
MOLS = os.environ.get("MOLS", "H2O OH").split()
DISSOC_DEFAULT = ("iface_N_D28:s1 iface_N_D28:s4 matrix_xp:s5 matrix_xp:s1 iface_N_D28:s2 matrix_xp:s4 iface_N_D28:s5")
os.makedirs("ads", exist_ok=True)


# ------------------------------------------------------------------ input handling
def load(path):
    s = open(path).read()
    head, rest = s.split("ATOMIC_SPECIES", 1)
    spec, rest = rest.split("CELL_PARAMETERS angstrom", 1)
    cell, rest = rest.split("ATOMIC_POSITIONS angstrom", 1)
    pos, tail = rest.split("K_POINTS", 1)
    cell = [list(map(float, l.split())) for l in cell.strip().splitlines()]
    atoms = []
    for l in pos.strip().splitlines():
        t = l.split()
        atoms.append((t[0], float(t[1]), float(t[2]), float(t[3])))
    species = [l.split() for l in spec.strip().splitlines()]
    return dict(head=head, species=species, cell=cell, atoms=atoms, tail="K_POINTS" + tail)


def mic_xy(d, L):
    return d - L * np.round(d / L)


def dist_to(p, atoms, cell):
    Lx, Ly = cell[0][0], cell[1][1]
    P = np.array([[a[1], a[2], a[3]] for a in atoms])
    d = P - p
    d[:, 0] = mic_xy(d[:, 0], Lx)
    d[:, 1] = mic_xy(d[:, 1], Ly)
    return np.sqrt((d ** 2).sum(1))


def write_input(fn, tmpl, atoms, free, prefix, calc="relax", ecut=None, extra_bands=0,
                new_species=(), pseudo_dir=None, extra_control="", gas_box=None):
    head = tmpl["head"]
    head = re.sub(r"calculation = '\w+'", f"calculation = '{calc}'", head)
    head = re.sub(r"prefix = '[^']*'", f"prefix = '{prefix}'", head)
    head = re.sub(r"outdir = '[^']*'", f"outdir = './tmp_{prefix}'", head)
    head = re.sub(r"tprnfor = \S+", "tprnfor = .true.", head)
    if calc == "relax":
        head = head.replace("&control\n", "&control\n  nstep = 80\n  forc_conv_thr = 2.0d-3\n  etot_conv_thr = 2.0d-4\n", 1)
    head = re.sub(r"nat = \d+", f"nat = {len(atoms)}", head)
    species = [list(s) for s in tmpl["species"]]
    have = {s[0] for s in species}
    for el in new_species:
        if el not in have:
            species.append([el, str(MASS[el]), PSEUDO_FILES[el]])
    head = re.sub(r"ntyp = \d+", f"ntyp = {len(species)}", head)
    if ecut:
        head = re.sub(r"ecutwfc = [\d.]+", f"ecutwfc = {ecut[0]}", head)
        head = re.sub(r"ecutrho = [\d.]+", f"ecutrho = {ecut[1]}", head)
    if extra_bands:
        head = re.sub(r"nbnd = (\d+)", lambda m: f"nbnd = {int(m.group(1)) + extra_bands}", head)
    if pseudo_dir:
        head = re.sub(r"pseudo_dir = '[^']*'", f"pseudo_dir = '{pseudo_dir}'", head)
    if calc == "relax" and "&ions" not in head:
        head = head.replace("&electrons", "&electrons", 1)
        # &ions goes after the &electrons block
        i = head.index("&electrons")
        j = head.index("\n/\n", i) + 3
        head = head[:j] + "&ions\n  ion_dynamics = 'bfgs'\n/\n" + head[j:]
    out = head + "ATOMIC_SPECIES\n" + "".join(f"  {s[0]}  {s[1]}  {s[2]}\n" for s in species)
    out += "CELL_PARAMETERS angstrom\n" + "".join("  %.6f  %.6f  %.6f\n" % tuple(v) for v in tmpl["cell"])
    out += "ATOMIC_POSITIONS angstrom\n"
    for a, f in zip(atoms, free):
        flag = "  1 1 1" if f else "  0 0 0"
        out += "  %s  %.6f  %.6f  %.6f%s\n" % (a[0], a[1], a[2], a[3], flag if calc == "relax" else "")
    out += "\n" + tmpl["tail"]
    open(fn, "w").write(out)


# ------------------------------------------------------------------ sites
def pick_amorphous_sites(cell_name, tmpl):
    rows = list(csv.DictReader(open(f"analysis_{cell_name}.csv")))
    atoms = tmpl["atoms"]
    cand = []
    for r in rows:
        i = int(r["atom"]) - 1
        if int(r["bottom"]) == 0 and float(r["depth"]) < DEPTH_MAX and r["e_p_win"] not in ("", "nan"):
            assert atoms[i][0] == r["el"], "atom order mismatch between CSV and QE input"
            cand.append((float(r["e_p_win"]), i, r))
    cand.sort()
    Lx, Ly = tmpl["cell"][0][0], tmpl["cell"][1][1]
    picks = []
    if SITE_MODE == "element":
        rng = np.random.default_rng(7)
        byel = {}
        for e, i, r in cand:
            byel.setdefault(r["el"], []).append((e, i, r))
        for k in byel:
            byel[k] = [byel[k][j] for j in rng.permutation(len(byel[k]))]
        progress = True
        while len(picks) < N_SITES and progress:
            progress = False
            for el in ("Ca", "Zn", "Mg"):
                if len(picks) >= N_SITES:
                    break
                lst = byel.get(el, [])
                while lst:
                    e, i, r = lst.pop(0)
                    if all(np.hypot(mic_xy(atoms[i][1] - atoms[j][1], Lx), mic_xy(atoms[i][2] - atoms[j][2], Ly)) >= MIN_SEP
                           for _, j, _ in picks):
                        picks.append((e, i, r)); progress = True
                        break
        return [dict(site=f"s{n+1}", atom=i, el=r["el"], n_ZnCa=r["n_ZnCa"], q=r["q"], e_p_win=r["e_p_win"],
                     e_p_occ=r["e_p_occ"], depth=r["depth"], mode="atop") for n, (e, i, r) in enumerate(picks)]
    targets = np.linspace(0, len(cand) - 1, N_SITES)
    for t in targets:
        order = sorted(range(len(cand)), key=lambda k: abs(k - t))
        for k in order:
            e, i, r = cand[k]
            if i in [p[1] for p in picks]:
                continue
            ok = all(np.hypot(mic_xy(atoms[i][1] - atoms[j][1], Lx), mic_xy(atoms[i][2] - atoms[j][2], Ly)) >= MIN_SEP
                     for _, j, _ in picks)
            if ok:
                picks.append((e, i, r)); break
    return [dict(site=f"s{n+1}", atom=i, el=r["el"], n_ZnCa=r["n_ZnCa"], q=r["q"], e_p_win=r["e_p_win"],
                 e_p_occ=r["e_p_occ"], depth=r["depth"], mode="atop") for n, (e, i, r) in enumerate(picks)]


def pick_crystal_sites(tmpl):
    atoms = tmpl["atoms"]; cell = tmpl["cell"]
    P = np.array([[a[1], a[2], a[3]] for a in atoms])
    ztop = P[:, 2].max()
    top = [i for i in range(len(atoms)) if P[i, 2] > ztop - 1.0]
    sec = [i for i in range(len(atoms)) if ztop - 3.4 < P[i, 2] < ztop - 1.8]
    Lx, Ly = cell[0][0], cell[1][1]

    def rel(i, j):
        d = P[j] - P[i]; d[0] = mic_xy(d[0], Lx); d[1] = mic_xy(d[1], Ly); return d
    i0 = max(top, key=lambda i: (P[i, 0], P[i, 1]))
    i0 = [i for i in top if abs(P[i, 0] - Lx / 2) < 4 and abs(P[i, 1] - Ly / 2) < 4][0]   # away from the cell edge
    nb = [j for j in top if j != i0 and np.hypot(*rel(i0, j)[:2]) < 3.5]
    sites = [dict(site="atop", atom=i0, el=atoms[i0][0], mode="atop", n_ZnCa="0", q="", e_p_win="", e_p_occ="", depth="0")]
    found = {}
    for a in range(len(nb)):
        for b in range(a + 1, len(nb)):
            if np.hypot(*rel(nb[a], nb[b])[:2]) < 3.5:
                c = P[i0][:2] + (rel(i0, nb[a])[:2] + rel(i0, nb[b])[:2]) / 3.0
                below = [k for k in sec if np.hypot(mic_xy(P[k, 0] - c[0], Lx), mic_xy(P[k, 1] - c[1], Ly)) < 0.7]
                kind = "hcp" if below else "fcc"
                if kind not in found:
                    found[kind] = c
    for kind in ("hcp", "fcc"):
        if kind in found:
            sites.append(dict(site=kind, atom=i0, el="Mg", mode=kind, xy=found[kind], n_ZnCa="0", q="",
                              e_p_win="", e_p_occ="", depth="0"))
    return sites


# ------------------------------------------------------------------ adsorbates
def place(tmpl, site, mol):
    atoms = tmpl["atoms"]; cell = tmpl["cell"]
    P = np.array([[a[1], a[2], a[3]] for a in atoms])
    ztop = P[:, 2].max()
    i = site["atom"]
    if site["mode"] == "atop":
        x, y, z = P[i]; h = 2.3 if mol == "H2O" else 1.95
    else:
        x, y = site["xy"]; z = ztop; h = 1.1 if mol == "OH" else 2.3
    p = np.array([x, y, z + h])
    for _ in range(20):                                  # lift until no metal atom is closer than 2.0 A
        if dist_to(p, atoms, cell).min() >= 2.0:
            break
        p[2] += 0.1
    if mol == "H2O":
        ads = [("O",) + tuple(p), ("H", p[0] + 0.767, p[1], p[2] + 0.594), ("H", p[0] - 0.767, p[1], p[2] + 0.594)]
    else:
        ads = [("O",) + tuple(p), ("H", p[0], p[1], p[2] + 0.97)]
    d = dist_to(p, atoms, cell)
    free = [(d[k] <= FREE_RADIUS and P[k, 2] > P[:, 2].min() + 4.0) for k in range(len(atoms))] + [True] * len(ads)
    return atoms + ads, free


def gas_template(ref):
    t = load(ref)
    t["head"] = re.sub(r"occupations = 'smearing'\n\s*smearing = '\w+'\n\s*degauss = [\d.eE+-]+\n",
                       "occupations = 'fixed'\n", t["head"])
    t["head"] = re.sub(r"\s*nbnd = \d+\n", "\n", t["head"])
    t["head"] = re.sub(r"\s*(tefield|dipfield) = \S+\n|\s*edir = \d\n|\s*emaxpos = \S+\n|\s*eopreg = \S+\n|\s*eamp = \S+\n", "\n", t["head"])
    t["head"] = re.sub(r"\n\s*\n", "\n", t["head"])
    L = 14.0
    t["cell"] = [[L, 0, 0], [0, L, 0], [0, 0, L]]
    t["species"] = [s for s in t["species"] if s[0] in ("Mg",)][:0]
    t["head"] = re.sub(r"mixing_beta = [\d.]+", "mixing_beta = 0.4", t["head"])
    t["head"] = t["head"].replace("  mixing_mode = 'local-TF'\n", "")
    t["tail"] = "K_POINTS gamma\n"
    return t


def gas_atoms(mol):
    c = 7.0
    if mol == "H2O":
        return [("O", c, c, c), ("H", c + 0.767, c, c + 0.594), ("H", c - 0.767, c, c + 0.594)]
    return [("H", c - 0.375, c, c), ("H", c + 0.375, c, c)]


# ------------------------------------------------------------------ build
def build():
    index = {"gas": [], "crystal": [], "matrix_H2O": [], "matrix_OH": [], "hi": []}
    first = None
    for cell in CELLS:
        ref = f"qe_inputs/{cell}_wf.in"
        tmpl = load(ref)
        pdir = re.search(r"pseudo_dir = '([^']*)'", tmpl["head"]).group(1)
        for el, fn in PSEUDO_FILES.items():
            if not os.path.exists(os.path.join(pdir, fn)):
                cand = sorted(glob.glob(os.path.join(pdir, el + "[._]*")) + glob.glob(os.path.join(pdir, el.lower() + "[._]*")))
                print(f"WARNING: {fn} not in {pdir}; candidates: {[os.path.basename(c) for c in cand]}")
        sites = pick_crystal_sites(tmpl) if cell.startswith("crystal") else pick_amorphous_sites(cell, tmpl)
        with open(f"ads/sites_{cell}.csv", "w", newline="") as f:
            w = csv.writer(f); w.writerow(["site", "atom", "el", "mode", "n_ZnCa", "q", "e_p_win", "e_p_occ", "depth"])
            for s in sites:
                w.writerow([s["site"], s["atom"] + 1, s["el"], s["mode"], s["n_ZnCa"], s["q"], s["e_p_win"], s["e_p_occ"], s["depth"]])
        print(f"{cell}: {len(sites)} sites -> " + ", ".join(f"{s['site']}({s['el']}#{s['atom']+1})" for s in sites))
        for s in sites:
            for mol in MOLS:
                if cell.startswith("crystal") and s["mode"] != "atop" and mol == "H2O":
                    continue                      # H2O does not sit in a hollow; OH tests hcp vs fcc
                atoms, free = place(tmpl, s, mol)
                name = f"ads_{cell}_{s['site']}_{mol}"
                write_input(f"qe_inputs/{name}.in", tmpl, atoms, free, name, "relax", extra_bands=EXTRA_BANDS,
                            new_species=("H", "O"))
                grp = "crystal" if cell.startswith("crystal") else f"matrix_{mol}"
                index[grp].append(name)
                if cell.startswith("crystal") and s["site"] == "atop" and mol == "H2O":
                    hname = name + "_hi"
                    write_input(f"qe_inputs/{hname}.in", tmpl, atoms, free, hname, "relax", ecut=HI,
                                extra_bands=EXTRA_BANDS, new_species=("H", "O"))
                    write_input(f"qe_inputs/{cell}_clean_hi.in", tmpl, tmpl["atoms"], [False] * len(tmpl["atoms"]),
                                f"{cell}_clean_hi", "scf", ecut=HI)
                    index["hi"] += [f"{cell}_clean_hi", hname]
        first = first or tmpl
    gt = gas_template(f"qe_inputs/{CELLS[0]}_wf.in")
    for mol in ("H2O", "H2"):
        at = gas_atoms(mol)
        write_input(f"qe_inputs/gas_{mol}.in", gt, at, [True] * len(at), f"gas_{mol}", "relax", new_species=("H", "O"))
        write_input(f"qe_inputs/gas_{mol}_hi.in", gt, at, [True] * len(at), f"gas_{mol}_hi", "relax", ecut=HI,
                    new_species=("H", "O"))
        index["gas"] += [f"gas_{mol}", f"gas_{mol}_hi"]
    for g, names in index.items():
        open(f"ads/list_{g}.txt", "w").write("\n".join(names) + "\n")
        print(f"list_{g}.txt: {len(names)} jobs")


# ------------------------------------------------------------------ analyze
def last_energy(fn):
    e = None
    for l in open(fn, errors="ignore"):
        if l.startswith("!"):
            e = float(l.split("=")[1].split()[0])
    return e


def status(fn):
    t = open(fn, errors="ignore").read()
    steps = t.count("number of scf cycles")
    ok = ("bfgs converged" in t) or ("End of BFGS Geometry Optimization" in t)
    done = "JOB DONE" in t
    return steps, ok, done


def final_geometry(out_fn, in_fn):
    """Final O position: height above the highest slab atom and number of metal neighbours within 2.5 A (xy periodic)."""
    try:
        cell = load(in_fn)["cell"]
        Lx, Ly = cell[0][0], cell[1][1]
        t = open(out_fn, errors="ignore").read()
        blk = t.split("ATOMIC_POSITIONS")[-1].splitlines()[1:]
        at = []
        for l in blk:
            p = l.split()
            if len(p) < 4:
                break
            try:
                at.append((p[0], float(p[1]), float(p[2]), float(p[3])))
            except ValueError:
                break
        o = [a for a in at if a[0] == "O"]
        metal = [a for a in at if a[0] in ("Mg", "Zn", "Ca")]
        if not o or not metal:
            return float("nan"), -1, -1
        o = o[0]
        ztop = max(a[3] for a in metal)
        cn = 0
        for a in metal:
            dx, dy = mic_xy(a[1] - o[1], Lx), mic_xy(a[2] - o[2], Ly)
            if (dx * dx + dy * dy + (a[3] - o[3]) ** 2) ** 0.5 < 2.5:
                cn += 1
        nh = 0
        for a in at:
            if a[0] == "H":
                dx, dy = mic_xy(a[1] - o[1], Lx), mic_xy(a[2] - o[2], Ly)
                if (dx * dx + dy * dy + (a[3] - o[3]) ** 2) ** 0.5 < 1.25:
                    nh += 1
        return o[3] - ztop, cn, nh
    except Exception:
        return float("nan"), -1, -1


def analyze():
    def E(name):
        fn = f"qe_inputs/{name}.out"
        return last_energy(fn) if os.path.exists(fn) else None
    rows = []
    ALLCELLS = sorted(os.path.basename(f)[6:-4] for f in glob.glob("ads/sites_*.csv")) or CELLS
    cmap = {}
    if os.path.exists("ads/cref_map.csv"):
        cmap = {r["ads"]: r["cref"] for r in csv.DictReader(open("ads/cref_map.csv"))}
    for tag, suf in (("", ""), ("_hi", "_hi")):
        eh2o, eh2 = E(f"gas_H2O{suf}"), E(f"gas_H2{suf}")
        if eh2o is None or eh2 is None:
            print(f"gas references{suf or ''} missing"); continue
        for cell in ALLCELLS:
            eclean = E(f"{cell}_clean_hi") if suf else E(f"{cell}_wf")
            if eclean is None:
                continue
            sfn = f"ads/sites_{cell}.csv"
            sites = {r["site"]: r for r in csv.DictReader(open(sfn))} if os.path.exists(sfn) else {}
            for out in sorted(glob.glob(f"qe_inputs/ads_{cell}_*{suf}.out")):
                name = os.path.basename(out)[:-4]
                if suf == "" and name.endswith("_hi"):
                    continue
                m = re.match(rf"ads_{cell}_(\w+?)_(H2O|OH)(_hi)?$", name)
                if not m:
                    continue
                site, mol = m.group(1), m.group(2)
                e = last_energy(out)
                steps, conv, done = status(out)
                if e is None:
                    continue
                ref = eh2o if mol == "H2O" else eh2o - 0.5 * eh2
                ea = (e - eclean - ref) * RY
                s = sites.get(site, {})
                Em = float("nan")
                if name in cmap:
                    ec = E(cmap[name])
                    if ec is not None and "JOB DONE" in open(f"qe_inputs/{cmap[name]}.out", errors="ignore").read():
                        Em = (e - ec - ref) * RY
                zO, cn, nh = final_geometry(out, f"qe_inputs/{name}.in")
                try:
                    dz = "%+.3f" % (float(s.get("q", "")) - ZVAL[s.get("el", "")])
                except (ValueError, KeyError):
                    dz = ""
                rows.append(dict(zO=zO, cn=cn, nh=nh, dz=dz, Em=Em, cell=cell, site=site, mol=mol, cut="hi" if suf else "std", Eads=ea, conv=conv, steps=steps,
                                 el=s.get("el", ""), n=s.get("n_ZnCa", ""), q=s.get("q", ""), pw=s.get("e_p_win", ""),
                                 po=s.get("e_p_occ", "")))
    if not rows:
        print("no finished adsorption runs yet"); return
    with open("ads/ads_table.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["cell", "site", "mol", "cut", "Eads", "Em", "el", "dZ", "n_ZnCa", "zO", "cn", "nh"])
        for r in rows:
            w.writerow([r["cell"], r["site"], r["mol"], r["cut"], r["Eads"], r["Em"], r["el"], r["dz"], r["n"], r["zO"], r["cn"], r["nh"]])
    print(f"{'cell':14s} {'site':5s} {'mol':4s} {'cut':3s} {'E_ads/eV':>9s} E_matched conv steps  el  n_ZnCa  dZ      e_p_win e_p_occ  zO/A  CN(O) nH")
    for r in rows:
        print(f"{r['cell']:14s} {r['site']:5s} {r['mol']:4s} {r['cut']:3s} {r['Eads']:9.3f} {r['Em']:9.3f} {str(r['conv']):5s} {r['steps']:4d}  {r['el']:3s} {r['n']:>5s}  {r['dz']:7s} {r['pw'][:7]:7s} {r['po'][:7]:7s} {r['zO']:5.2f} {r['cn']:3d} {r['nh']:3d}")
    # ---- region summary: crystal / glass matrix / interphase (site-matched energies)
    def region(c):
        return "crystal" if c.startswith("crystal") else ("interphase" if c.startswith("iface") else "glass matrix")
    print("\nRegion summary of E_ads(site-matched) in eV (n sites; mean +- sd; min / median / max)")
    groups = {}
    for r in rows:
        if r["cut"] == "std" and not np.isnan(r["Em"]):
            groups.setdefault((region(r["cell"]), r["mol"]), []).append(r["Em"])
    for (reg, mol), v in sorted(groups.items()):
        v = np.array(v)
        sd = v.std(ddof=1) if len(v) > 1 else float("nan")
        print(f"  {reg:13s} {mol:3s} n={len(v):2d}  {v.mean():+.3f} +- {sd:.3f}   {v.min():+.3f} / {np.median(v):+.3f} / {v.max():+.3f}")
    rng = np.random.default_rng(1)

    def perm_p(x, y, n=20000):
        """two-sided permutation p-value for the difference of means (numpy only, no scipy)"""
        x, y = np.array(x), np.array(y)
        obs = abs(y.mean() - x.mean()); pool = np.concatenate([x, y]); k = len(x); c = 0
        for _ in range(n):
            q = rng.permutation(pool)
            c += abs(q[k:].mean() - q[:k].mean()) >= obs - 1e-12
        return (c + 1) / (n + 1)
    for mol in ("H2O", "OH"):
        ga, gb = groups.get(("glass matrix", mol)), groups.get(("interphase", mol))
        if ga and gb and len(ga) >= 3 and len(gb) >= 3:
            print(f"  {mol}: interphase minus glass matrix  mean {np.mean(gb) - np.mean(ga):+.3f} eV, "
                  f"sd {np.std(gb, ddof=1):.3f} vs {np.std(ga, ddof=1):.3f}, permutation p = {perm_p(ga, gb):.3f}  (n = {len(gb)} vs {len(ga)})")
    print("\nBy the metal atom under the adsorbate (amorphous cells, site-matched, eV; mean +- sd, n)")
    byel = {}
    for r in rows:
        if r["cut"] == "std" and not np.isnan(r["Em"]) and not r["cell"].startswith("crystal"):
            byel.setdefault((r["mol"], r["el"]), []).append(r["Em"])
    for (mol, el), v in sorted(byel.items()):
        sd = np.std(v, ddof=1) if len(v) > 1 else float("nan")
        print(f"  {mol:3s} {el}:  {np.mean(v):+.3f} +- {sd:.3f}   n={len(v)}")
    nb = [f"{r['cell']}:{r['site']}" for r in rows if r["cut"] == "std" and r["mol"] == "H2O" and r["cn"] == 0]
    if nb:
        print("  H2O with no metal neighbour within 2.5 A (physisorbed, not bonded):", ", ".join(nb))
    print()
    # ---- water dissociation (H2O* -> OH* + H*), same free atoms in both states
    drows = []
    for fn in sorted(glob.glob("qe_inputs/ddis_*.out")):
        base = os.path.basename(fn)[5:-4]
        ref = f"qe_inputs/dref_{base}.out"
        if not os.path.exists(ref) or "JOB DONE" not in open(fn, errors="ignore").read() \
                or "JOB DONE" not in open(ref, errors="ignore").read():
            continue
        site = base.split("_")[-1]; cell = base[:-(len(site) + 1)]
        zO, cn, nh = final_geometry(fn, f"qe_inputs/ddis_{base}.in")
        ea = [r["Em"] for r in rows if (r["cell"], r["site"], r["mol"], r["cut"]) == (cell, site, "H2O", "std")]
        el = next((r["el"] for r in rows if (r["cell"], r["site"]) == (cell, site)), "")
        e1 = ea[0] if ea else float("nan")
        de = (last_energy(fn) - last_energy(ref)) * RY
        drows.append((cell, site, el, e1, de, e1 + de, "recombined to H2O" if nh == 2 else "dissociated"))
    if drows:
        print("Water dissociation H2O* -> OH* + H* (eV; negative = downhill; thermodynamics only, no barrier)")
        print("  E_final = E_ads(H2O) + dE_diss = energy of OH* + H* relative to gas H2O (mixes two references; indicative)")
        print(f"  {'cell':14s} {'site':4s} {'el':2s} {'E_ads(H2O)':>10s} {'dE_diss':>8s} {'E_final':>8s}  final state")
        for c, st, el, e1, e2, ef, flag in drows:
            print(f"  {c:14s} {st:4s} {el:2s} {e1:10.3f} {e2:8.3f} {ef:8.3f}  {flag}")
        good = [d for d in drows if d[6] == "dissociated"]
        for el in ("Ca", "Mg", "Zn"):
            v = [d[5] for d in good if d[2] == el]
            if v:
                print(f"  E_final {el}: mean {np.mean(v):+.3f}  n={len(v)}")
        print()
    # cutoff test
    for r in rows:
        if r["cut"] == "hi":
            lo = [x for x in rows if x["cut"] == "std" and (x["cell"], x["site"], x["mol"]) == (r["cell"], r["site"], r["mol"])]
            if lo:
                print(f"cutoff test {r['cell']} {r['site']} {r['mol']}: {lo[0]['Eads']:.3f} eV at 40/200 vs {r['Eads']:.3f} eV at 60/480 "
                      f"(difference {lo[0]['Eads'] - r['Eads']:+.3f} eV)")
    # descriptor correlations on the amorphous cells
    def rank(a): return np.argsort(np.argsort(a)).astype(float)
    for ekey, elab in (("Eads", "E_ads"), ("Em", "E_ads(site-matched)")):
        for mol in ("H2O", "OH"):
            for desc, key in (("e_p_win", "pw"), ("e_p_occ", "po"), ("dZ", "dz"), ("n_ZnCa", "n")):
                sel = [r for r in rows if r["cut"] == "std" and r["mol"] == mol and r[key] not in ("", None)
                       and not r["cell"].startswith("crystal") and not np.isnan(r[ekey])]
                if len(sel) >= 4:
                    x = np.array([float(r[key]) for r in sel]); y = np.array([r[ekey] for r in sel])
                    rx, ry = rank(x), rank(y)
                    rho = np.corrcoef(rx, ry)[0, 1]
                    cnt = sum(abs(np.corrcoef(rx, rng.permutation(ry))[0, 1]) >= abs(rho) - 1e-12 for _ in range(5000))
                    print(f"{mol}: Spearman rho({elab}, {desc}) = {rho:+.2f}  (n={len(sel)}, permutation p = {(cnt + 1) / 5001:.3f})")


# ------------------------------------------------------------------ water dissociation energy
def final_positions(out_fn):
    t = open(out_fn, errors="ignore").read()
    blk = t.split("ATOMIC_POSITIONS")[-1].splitlines()[1:]
    at = []
    for l in blk:
        p = l.split()
        if len(p) < 4:
            break
        try:
            at.append((p[0], float(p[1]), float(p[2]), float(p[3])))
        except ValueError:
            break
    return at


def dissoc():
    """H2O* -> OH* + H* on finished H2O adsorption sites (DISSOC_SITES="cell:site ..."). For each site two relaxations
    with the SAME free atoms: dref_<cell>_<site> (molecular water, restarted from the relaxed structure) and
    ddis_<cell>_<site> (one H moved to the nearest surface metal atom)."""
    lists = [[], [], []]
    k = 0
    for item in os.environ.get("DISSOC_SITES", DISSOC_DEFAULT).split():
        cell, site = item.split(":")
        name = f"ads_{cell}_{site}_H2O"
        out = f"qe_inputs/{name}.out"
        if not (os.path.exists(out) and "JOB DONE" in open(out, errors="ignore").read()):
            print(f"skip {item}: {name} not finished"); continue
        tmpl = load(f"qe_inputs/{cell}_wf.in"); n = len(tmpl["atoms"]); cp = tmpl["cell"]
        Lx, Ly = cp[0][0], cp[1][1]
        fin = final_positions(out)
        if len(fin) != n + 3:
            print(f"skip {item}: unexpected atom count in {out}"); continue
        slab, ads = fin[:n], fin[n:]
        free_old = [a[4] for a in parse_flags(f"qe_inputs/{name}.in")][:n]
        O = np.array(ads[0][1:]); Hs = [np.array(a[1:]) for a in ads[1:]]
        mv = 0 if Hs[0][2] <= Hs[1][2] else 1                  # the H closer to the surface moves
        P = np.array([[a[1], a[2], a[3]] for a in slab])

        def lat(p, q):
            return float(np.hypot(mic_xy(q[0] - p[0], Lx), mic_xy(q[1] - p[1], Ly)))
        near = [j for j in range(n) if lat(O, P[j]) < 3.4]
        ztl = max(P[j, 2] for j in near) if near else P[:, 2].max()
        cands = [j for j in near if P[j, 2] > ztl - 1.2 and lat(O, P[j]) > 1.8]
        if cands:
            kk = min(cands, key=lambda j: lat(Hs[mv], P[j]))
            dx, dy = mic_xy(P[kk, 0] - O[0], Lx), mic_xy(P[kk, 1] - O[1], Ly)
            tgt = np.array([O[0] + dx * 0.85, O[1] + dy * 0.85, P[kk, 2] + 1.8])
        else:
            tgt = np.array([O[0] + 2.4, O[1], ztl + 1.0])
        for _ in range(25):
            if dist_to(tgt, slab, cp).min() >= 1.6:
                break
            tgt[2] += 0.1
        newads = [ads[0], ads[2 - mv], ("H", float(tgt[0]), float(tgt[1]), float(tgt[2]))]
        d = dist_to(tgt, slab, cp); zmin = P[:, 2].min()
        free2 = [bool(free_old[j] or (d[j] <= 4.0 and P[j, 2] > zmin + 4.0)) for j in range(n)] + [True] * 3
        base = f"{cell}_{site}"
        write_input(f"qe_inputs/dref_{base}.in", tmpl, slab + ads, free2, f"dref_{base}", "relax",
                    extra_bands=EXTRA_BANDS, new_species=("H", "O"))
        write_input(f"qe_inputs/ddis_{base}.in", tmpl, slab + newads, free2, f"ddis_{base}", "relax",
                    extra_bands=EXTRA_BANDS, new_species=("H", "O"))
        lists[k % 3] += [f"dref_{base}", f"ddis_{base}"]; k += 1
        print(f"{item}: H moved to ({tgt[0]:.2f},{tgt[1]:.2f},{tgt[2]:.2f}), {sum(free2)} free atoms")
    for q in range(3):
        open(f"ads/list_diss_{q+1}.txt", "w").write("\n".join(lists[q]) + ("\n" if lists[q] else ""))
    print("dissociation lists:", [len(l) for l in lists])


# ------------------------------------------------------------------ site-matched clean references
def parse_flags(fn):
    t = open(fn).read()
    pos = t.split("ATOMIC_POSITIONS angstrom", 1)[1].split("K_POINTS", 1)[0]
    out = []
    for l in pos.strip().splitlines():
        p = l.split()
        out.append((p[0], float(p[1]), float(p[2]), float(p[3]), len(p) >= 7 and p[4] == "1"))
    return out


def cleanref():
    """For every ads_<cell>_<site>_<mol>.in: the same slab with the same free atoms but WITHOUT adsorbate, relaxed.
    E_ads(site-matched) = E(slab+ads) - E(slab, same free atoms, relaxed) - E(gas reference).
    Entries of cells that are not in CELLS stay in ads/cref_map.csv."""
    old = {}
    if os.path.exists("ads/cref_map.csv"):
        old = {r["ads"]: r["cref"] for r in csv.DictReader(open("ads/cref_map.csv"))}
    new = {}
    for cell in CELLS:
        tmpl = load(f"qe_inputs/{cell}_wf.in")
        n = len(tmpl["atoms"])
        sets, mine = {}, []
        for fn in sorted(glob.glob(f"qe_inputs/ads_{cell}_*.in")):
            name = os.path.basename(fn)[:-3]
            if name.endswith("_hi"):
                continue
            at = parse_flags(fn)
            slab = at[:n]
            assert all(a[0] in ("Mg", "Zn", "Ca") for a in slab) and len(at) > n, f"unexpected atom order in {fn}"
            assert all(abs(a[1] - b[1]) < 1e-4 and abs(a[2] - b[2]) < 1e-4 and abs(a[3] - b[3]) < 1e-4
                       for a, b in zip(slab, tmpl["atoms"])), f"slab positions differ from {cell}_wf.in in {fn}"
            key = tuple(a[4] for a in slab)
            if key not in sets:
                cname = f"cref_{cell}_{len(sets) + 1}"
                sets[key] = cname
                write_input(f"qe_inputs/{cname}.in", tmpl, tmpl["atoms"], list(key), cname, "relax")
            new[name] = sets[key]
        names = list(sets.values())
        open(f"ads/list_cref_{cell}_a.txt", "w").write("\n".join(names[0::2]) + "\n")
        open(f"ads/list_cref_{cell}_b.txt", "w").write("\n".join(names[1::2]) + ("\n" if names[1::2] else ""))
        print(f"{cell}: {len(new) and sum(1 for k in new if k.startswith('ads_' + cell + '_'))} adsorption runs share {len(sets)} clean references")
    old.update(new)
    with open("ads/cref_map.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["ads", "cref"]); w.writerows(sorted(old.items()))
    # legacy lists (crystal / matrix groups) kept for the first calibration
    crys = sorted({c for a_, c in old.items() if c.startswith("cref_crystal")})
    mat = sorted({c for a_, c in old.items() if c.startswith("cref_matrix")})
    open("ads/list_cref_crystal.txt", "w").write("\n".join(crys) + "\n")
    for k in range(3):
        open(f"ads/list_cref_matrix_{k+1}.txt", "w").write("\n".join(mat[k::3]) + ("\n" if mat[k::3] else ""))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "build":
        build()
    elif cmd == "analyze":
        analyze()
    elif cmd == "cleanref":
        cleanref()
    elif cmd == "dissoc":
        dissoc()
    else:
        print(__doc__)
