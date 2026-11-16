"""
Tier-1 referee-proofing checks (builds inputs only; run them with submit_checks.sh).

  python3 16_checks.py build     # writes every check input listed below
  python3 16_checks.py analyze   # compares each check against the production value

  kpt   <cell>_k221        : the dipole-corrected clean slab of one amorphous cell at 2x2x1 instead of Gamma.
                             Checks the k-point convergence that was so far only verified on the crystal.
  cut   *_hi2              : one amorphous H2O site + its site-matched clean reference at 60/480 Ry.
                             The existing 60/480 test is crystal-only; gas_H2O_hi is reused.
  d3    *_d3               : three H2O sites (one per element) + their references + the gas, with
                             vdw_corr='grimme-d3'. Tests whether the Ca > Mg > Zn ordering survives dispersion,
                             which the manuscript currently assumes.
  wffix iface_N_D27_wf     : rerun of the one work function that never converged, with a gentler mixer.

Env: KCELL (default matrix_xp), CUTSITE (default matrix_xp:s3), D3SITES (default one site per element).
NOTE: 'grimme-d3' must be available in your pw.x build. If the d3 runs die in iosys, drop them and say so in
the manuscript instead of guessing the offset.
"""
import csv, importlib, os, re, sys

A = importlib.import_module("11_adsorption")
RY = A.RY
KCELL = os.environ.get("KCELL", "matrix_xp")
CUTSITE = os.environ.get("CUTSITE", "matrix_xp:s3")
D3SITES = os.environ.get("D3SITES", "matrix_xp:s1 iface_N_D28:s4 iface_N_D28:s5").split()
HI = (60.0, 480.0)


def cmap():
    f = "ads/cref_map.csv"
    return {r["ads"]: r["cref"] for r in csv.DictReader(open(f))} if os.path.exists(f) else {}


def retune(src, dst, prefix, subs):
    """copy a QE input, rename its prefix/outdir, and apply (pattern, replacement) pairs in &system/&electrons"""
    s = open(src).read()
    old = re.search(r"prefix = '([^']*)'", s).group(1)
    s = s.replace(f"prefix = '{old}'", f"prefix = '{prefix}'")
    s = re.sub(r"outdir = '[^']*'", f"outdir = './tmp_{prefix}'", s)
    for pat, rep in subs:
        s, n = re.subn(pat, rep, s, count=1)
        if n == 0 and rep:
            s = s.replace("&system\n", "&system\n  " + rep.strip() + "\n", 1)
    open(dst, "w").write(s)
    print(f"wrote {dst}")


def build():
    names = {"kpt": [], "cut": [], "d3": [], "wffix": []}
    # ---- k-points on one amorphous cell
    src = f"qe_inputs/{KCELL}_wf.in"
    if os.path.exists(src):
        n = f"{KCELL}_wf_k221"
        s = open(src).read()
        s = s.replace(f"prefix = '{KCELL}_wf'", f"prefix = '{n}'")
        s = re.sub(r"outdir = '[^']*'", f"outdir = './tmp_{n}'", s)
        s = re.sub(r"K_POINTS\s+\w+\s*\n[^\n]*\n?", "K_POINTS automatic\n  2 2 1 0 0 0\n", s)
        if "K_POINTS automatic" not in s:
            s = re.sub(r"K_POINTS.*", "K_POINTS automatic\n  2 2 1 0 0 0", s, flags=re.S)
        open(f"qe_inputs/{n}.in", "w").write(s)
        print(f"wrote qe_inputs/{n}.in"); names["kpt"].append(n)
    else:
        print(f"skip kpt: {src} missing")
    # ---- cutoff on one amorphous adsorption site
    cm = cmap()
    cell, site = CUTSITE.split(":")
    base = f"ads_{cell}_{site}_H2O"
    hi = [(r"ecutwfc = [\d.]+", f"ecutwfc = {HI[0]}"), (r"ecutrho = [\d.]+", f"ecutrho = {HI[1]}")]
    for nm in (base, cm.get(base)):
        if nm and os.path.exists(f"qe_inputs/{nm}.in"):
            retune(f"qe_inputs/{nm}.in", f"qe_inputs/{nm}_hi2.in", f"{nm}_hi2", hi)
            names["cut"].append(f"{nm}_hi2")
        elif nm:
            print(f"skip cut: qe_inputs/{nm}.in missing")
    # ---- dispersion on one site per element
    d3 = [(r"\n", "\n  vdw_corr = 'grimme-d3'\n", )]
    for item in D3SITES:
        cell, site = item.split(":")
        b = f"ads_{cell}_{site}_H2O"
        for nm in (b, cm.get(b)):
            if nm and os.path.exists(f"qe_inputs/{nm}.in"):
                retune(f"qe_inputs/{nm}.in", f"qe_inputs/{nm}_d3.in", f"{nm}_d3",
                       [(r"vdw_corr *= *'[^']*'", "vdw_corr = 'grimme-d3'")])
                names["d3"].append(f"{nm}_d3")
            elif nm:
                print(f"skip d3: qe_inputs/{nm}.in missing")
    if os.path.exists("qe_inputs/gas_H2O.in"):
        retune("qe_inputs/gas_H2O.in", "qe_inputs/gas_H2O_d3.in", "gas_H2O_d3",
               [(r"vdw_corr *= *'[^']*'", "vdw_corr = 'grimme-d3'")])
        names["d3"].append("gas_H2O_d3")
    # ---- the failed work function, gentler mixer
    src = "qe_inputs/iface_N_D27_wf.in"
    if os.path.exists(src):
        retune(src, "qe_inputs/iface_N_D27_wf_fix.in", "iface_N_D27_wf_fix",
               [(r"mixing_beta *= *[\d.]+", "mixing_beta = 0.05")])
        s = open("qe_inputs/iface_N_D27_wf_fix.in").read()
        s = re.sub(r"\n *electron_maxstep *= *\d+", "", s)
        s = s.replace("&electrons\n", "&electrons\n  electron_maxstep = 400\n  mixing_ndim = 12\n", 1)
        open("qe_inputs/iface_N_D27_wf_fix.in", "w").write(s)
        names["wffix"].append("iface_N_D27_wf_fix")
    for k, v in names.items():
        if v:
            open(f"ads/list_chk_{k}.txt", "w").write("\n".join(v) + "\n")
    print("\nlists:", {k: len(v) for k, v in names.items()})


def E(n):
    f = f"qe_inputs/{n}.out"
    return A.last_energy(f) if os.path.exists(f) else None


def analyze():
    cm = cmap()
    print("k-point check (clean dipole-corrected slab, Gamma vs 2x2x1)")
    a, b = E(f"{KCELL}_wf"), E(f"{KCELL}_wf_k221")
    print(f"  {KCELL}: {'%.4f eV/atom difference' % ((b - a) * RY) if a and b else 'pending'}")
    print("\ncutoff check on an amorphous site (40/200 vs 60/480 Ry)")
    cell, site = CUTSITE.split(":")
    base = f"ads_{cell}_{site}_H2O"
    for tag, suf in (("40/200", ""), ("60/480", "_hi2")):
        e, c = E(base + suf), E((cm.get(base) or "") + suf)
        g = E("gas_H2O_hi" if suf else "gas_H2O")
        print(f"  {tag}: " + (f"E_ads = {(e - c - g) * RY:+.3f} eV" if None not in (e, c, g) else "pending"))
    print("\ndispersion check (PBE vs PBE+D3), site-matched")
    for item in D3SITES:
        cl, st = item.split(":")
        b = f"ads_{cl}_{st}_H2O"
        row = []
        for tag, suf, gas in (("PBE", "", "gas_H2O"), ("PBE+D3", "_d3", "gas_H2O_d3")):
            e, c, g = E(b + suf), E((cm.get(b) or "") + suf), E(gas)
            row.append(f"{(e - c - g) * RY:+.3f}" if None not in (e, c, g) else "pending")
        print(f"  {item:22s} PBE {row[0]:>8s}   PBE+D3 {row[1]:>8s}")
    print("\nwork function rerun")
    print("  iface_N_D27: " + ("converged" if E("iface_N_D27_wf_fix") else "pending"))


if __name__ == "__main__":
    {"build": build, "analyze": analyze}.get(sys.argv[1] if len(sys.argv) > 1 else "", lambda: print(__doc__))()
