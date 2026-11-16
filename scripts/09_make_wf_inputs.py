"""Work-function inputs: copy qe_inputs/<name>.in -> qe_inputs/<name>_wf.in with a slab dipole correction.
Everything else (cutoffs, smearing, nbnd, Gamma point) is identical to the production SCF, so the
only difference is tefield/dipfield. Usage: python3 09_make_wf_inputs.py NAME [NAME ...]"""
import re, sys

def make(name):
    s = open(f"qe_inputs/{name}.in").read()
    # cell height and highest atom
    cp = re.search(r"CELL_PARAMETERS angstrom\s*\n(.*?)\n\s*\n?ATOMIC_POSITIONS", s, re.S)
    rows = [list(map(float, l.split())) for l in cp.group(1).strip().splitlines()]
    cz = rows[2][2]
    pos = s.split("ATOMIC_POSITIONS angstrom")[1].split("K_POINTS")[0].strip().splitlines()
    zs = [float(l.split()[3]) for l in pos if l.strip()]
    ztop, zbot = max(zs), min(zs)
    vac = cz - (ztop - zbot)
    emaxpos = (ztop + 0.5 * (cz - ztop + zbot)) / cz       # middle of the vacuum gap
    emaxpos = min(max(emaxpos, 0.05), 0.95)
    s = s.replace(f"prefix = '{name}'", f"prefix = '{name}_wf'").replace(f"./tmp_{name}'", f"./tmp_{name}_wf'")
    add = (f"  tefield = .true.\n  dipfield = .true.\n  edir = 3\n  emaxpos = {emaxpos:.3f}\n"
           f"  eopreg = 0.08\n  eamp = 0.0\n")
    s = s.replace("/\n&electrons", add + "/\n&electrons", 1)
    open(f"qe_inputs/{name}_wf.in", "w").write(s)
    print(f"{name}: cell z = {cz:.2f} A, slab {zbot:.2f}..{ztop:.2f} A, vacuum {vac:.1f} A, emaxpos = {emaxpos:.3f}")

for n in sys.argv[1:]:
    make(n)
