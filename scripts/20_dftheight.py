"""
Does DFT agree with MEAM that Zn sits higher than Mg at the free surface?  (no new runs, seconds)

19_relaxprofile.py: under MEAM every element relaxes INWARD (Mg -0.61, Ca -0.37, Zn -0.22 A), so Zn
ends 0.40 A above Mg and enriches the outermost layer. The sign is what this potential's known Zn
underbinding manufactures - weaker binding, less inward pull - so the magnitude needs an independent check.

It is already on disk. The cref_* clean slabs were DFT-relaxed from the MEAM geometry with the same free
atoms, so their outputs hold the initial and final coordinates. This parses both and measures, per atom,
how much further DFT moves each element:
  d(Zn) - d(Mg) > 0   DFT splits them further  -> MEAM understates it, the claim is safe
  d(Zn) - d(Mg) ~ 0   DFT is content           -> MEAM geometry is at a DFT stationary point, claim stands
  d(Zn) - d(Mg) < 0   DFT pulls Zn back in     -> MEAM exaggerated it; scale the claim, or drop the
                                                  0.041 eV offset that the composition comparison rests on
Only free atoms count: the fixed slab bottom cannot move and would dilute the means toward zero.

  python3 20_dftheight.py
  GLOB='qe_inputs/cref_*.out' python3 20_dftheight.py
"""
import glob, os, re, sys
import numpy as np

ELS = ("Ca", "Mg", "Zn")
rng = np.random.default_rng(0)
PAT = os.environ.get("GLOB", "qe_inputs/cref_*.out")


def blocks(fn):
    """every ATOMIC_POSITIONS block in a QE output, as (elements, z, free_mask)"""
    t = open(fn, errors="ignore").read()
    out = []
    for m in re.finditer(r"ATOMIC_POSITIONS[^\n]*\n", t):
        el, z, fr = [], [], []
        for l in t[m.end():].splitlines():
            p = l.split()
            if len(p) < 4 or p[0] not in ELS + ("H", "O"):
                break
            el.append(p[0]); z.append(float(p[3]))
            fr.append(not (len(p) >= 7 and p[4] == "0" and p[5] == "0" and p[6] == "0"))
        if el:
            out.append((np.array(el), np.array(z), np.array(fr)))
    return out


def perm(d):
    d = np.asarray(d, float); d = d[~np.isnan(d)]
    if len(d) < 3:
        return np.nan, np.nan, np.nan, 0
    m = d.mean()
    c = sum(abs((d * rng.choice([-1, 1], len(d))).mean()) >= abs(m) - 1e-12 for _ in range(20000))
    return m, d.std(ddof=1) / np.sqrt(len(d)), (c + 1) / 20001, len(d)


cells, skipped = [], []
for fn in sorted(glob.glob(PAT)):
    nm = os.path.basename(fn)[:-4]
    b = blocks(fn)
    if len(b) < 2:
        skipped.append((nm, f"{len(b)} position block(s) - not a relaxation, or still running")); continue
    (e0, z0, f0), (e1, z1, f1) = b[0], b[-1]
    if len(e0) != len(e1) or not (e0 == e1).all():
        skipped.append((nm, "atom list changed between first and last block")); continue
    if f0.sum() == 0:
        skipped.append((nm, "no free atoms")); continue
    # reference out the rigid shift of the whole free set, so this measures splitting, not translation
    d = (z1 - z0)
    cells.append((nm, e0, d, f0, len(b) - 1))
if skipped:
    print("skipped:")
    for nm, why in skipped:
        print(f"  {nm:34s} {why}")
    print()
if not cells:
    print(f"no relaxed cref outputs matched {PAT}"); sys.exit()

print(f"{len(cells)} DFT-relaxed clean slabs, free atoms only\n")
print(f"{'cell':34s} {'steps':>5s} {'free':>5s} | " + " ".join(f"{'d' + e:>7s}" for e in ELS) + f" {'dZn-dMg':>9s}")
per = {e: [] for e in ELS}
rel = []
for nm, el, d, fr, ns in cells:
    mu = {}
    for e in ELS:
        m = fr & (el == e)
        mu[e] = d[m].mean() if m.sum() else np.nan
        if m.sum():
            per[e].append(mu[e])
    r = mu["Zn"] - mu["Mg"]
    if not np.isnan(r):
        rel.append(r)
    print(f"{nm:34s} {ns:5d} {int(fr.sum()):5d} | "
          + " ".join(f"{mu[e]:+7.3f}" if not np.isnan(mu[e]) else f"{'-':>7s}" for e in ELS)
          + (f" {r:+9.3f}" if not np.isnan(r) else f" {'-':>9s}"))

print("\nmean DFT displacement of the free atoms (A, + = outward)")
for e in ELS:
    m, se, p, k = perm(per[e])
    print(f"  {e:>3s} {m:+6.3f} +- {se:5.3f}   p = {p:.4f}   ({k} cells)")
m, se, p, k = perm(rel)
print(f"\nZn relative to Mg under DFT: {m:+.3f} +- {se:.3f} A   paired p = {p:.4f}   ({k} cells)")
print(f"MEAM gave                    +0.400 +- 0.040 A")
if not np.isnan(m):
    if m > 2 * se and m > 0:
        v = "DFT splits Zn and Mg FURTHER -> MEAM understates the splitting; the claim is safe."
    elif abs(m) <= 2 * se:
        v = ("DFT leaves the splitting alone -> the MEAM geometry is at a DFT stationary point.\n"
             "  The 0.4 A splitting stands, but say it was checked against DFT and give this number.")
    else:
        v = ("DFT pulls Zn back IN -> MEAM exaggerated the splitting.\n"
             "  Report the enrichment as a model-dependent upper bound, and do not rest the composition\n"
             "  comparison on the 0.041 eV offset without recomputing it from a DFT-relaxed surface.")
    print("\n  " + v)
print("""
Caveat to state either way: these relaxations freed only the top layers for a short BFGS, so the
absolute displacements are small and it is the Zn-vs-Mg difference, not the magnitude, that is tested.""")
