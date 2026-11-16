"""Work function from a dipole-corrected slab SCF:  W = V_vacuum(above the top surface) - E_F.
Needs  wf/<name>.cube  (pp.x plot_num=11, V_bare+V_H incl. the dipole sawtooth)  and
qe_inputs/<name>_wf.out.  Usage: python3 10_workfunction.py NAME [NAME ...]
pp.x writes the potential in Ry for plot_num=11, but the cube convention is Hartree, so BOTH
readings are printed; the one that reproduces Mg(0001) (~3.6-3.7 eV) in crystal_ideal is the right one."""
import re, sys
import numpy as np

BOHR = 0.529177210903
RY, HA = 13.605693123, 27.211386246

def read_cube(fn):
    with open(fn) as f:
        f.readline(); f.readline()
        nat = int(f.readline().split()[0])
        n, v = [], []
        for _ in range(3):
            t = f.readline().split(); n.append(int(t[0])); v.append([float(x) for x in t[1:4]])
        for _ in range(abs(nat)):
            f.readline()
        data = np.fromstring(f.read(), sep=" ")
    return np.array(n), np.array(v), data.reshape(n[0], n[1], n[2])

def fermi(fn):
    ef = None
    for l in open(fn, errors="ignore"):
        if "the Fermi energy is" in l:
            ef = float(l.split()[4])
    return ef

def top_z(name):
    s = open(f"qe_inputs/{name}_wf.in").read()
    pos = s.split("ATOMIC_POSITIONS angstrom")[1].split("K_POINTS")[0].strip().splitlines()
    return max(float(l.split()[3]) for l in pos if l.strip())

def run(name):
    ef = fermi(f"qe_inputs/{name}_wf.out")
    n, v, V = read_cube(f"wf/{name}.cube")
    prof = V.mean(axis=(0, 1))
    dz = v[2][2] * BOHR
    z = np.arange(n[2]) * dz
    lz = n[2] * dz
    zt = top_z(name)
    lo, hi = zt + 4.0, min(zt + 9.0, lz - 3.0)
    m = (z >= lo) & (z <= hi)
    if m.sum() < 5:
        print(name, "ERROR: vacuum window too thin"); return
    p = prof[m]
    slope = np.polyfit(z[m], p, 1)[0]
    print(f"\n{name}: E_F = {ef:.4f} eV | top atom z = {zt:.2f} A | vacuum window {lo:.1f}-{hi:.1f} A of {lz:.1f} A")
    for label, f in [("potential read as Ry", RY), ("potential read as Hartree", HA)]:
        vv = p.mean() * f
        print(f"  {label:27s}: V_vac = {vv:8.3f} eV (sd {p.std()*f:.3f}, slope {slope*f:+.4f} eV/A)  ->  W = {vv - ef:6.3f} eV")
    # step between the two vacuum sides as a sanity check on the dipole correction
    zb = np.where(z > zt + 1.0)[0]
    print("  profile every 2 A above the slab (Ry units):", " ".join(f"{prof[np.argmin(abs(z-(zt+d)))]:.4f}" for d in range(1, int(min(12, lz-zt-1)), 2)))

for nme in sys.argv[1:]:
    run(nme)
