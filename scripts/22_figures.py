"""
Publication figures and tables for the merged surface paper (22_figures.py).

  python3 22_figures.py inspect   # dump every source schema, check matplotlib. Writes nothing.
  python3 22_figures.py make      # build what the data supports -> figs/*.pdf, figs/*.png, figs/table*.csv

Each figure is independent: one that lacks its source is skipped with a reason, the rest still build.
A closing summary says what was made and what was missing.

matplotlib in apps/anaconda3/2024.06 is compiled against numpy 1.x, while ~/.local holds numpy 2.2.6 and
shadows it, which is the same clash that makes scipy unusable. Run with PYTHONNOUSERSITE=1 (submit_figures.sh
sets it) so ~/.local is ignored and anaconda's own numpy is used. Nothing here needs numpy 2.
"""
import csv, glob, os, re, sys
import numpy as np

EL = ("Ca", "Mg", "Zn")
C = {"Ca": "#D55E00", "Mg": "#0072B2", "Zn": "#009E73"}
MK = {"interphase": "s", "glass matrix": "o", "crystal": "^", "series glass": "o"}
made, skipped = [], []


def region(cell):
    if cell.startswith("iface_"):
        return "interphase"
    if cell.startswith("matrix_"):
        return "glass matrix"
    if cell.startswith("crystal"):
        return "crystal"
    return "series glass"


def rows(fn):
    return list(csv.DictReader(open(fn))) if os.path.exists(fn) else []


def F(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def newest(pat):
    f = sorted(glob.glob(pat), key=os.path.getmtime)
    return f[-1] if f else None


def ads(mol=None, amorphous_only=True, cut="std"):
    out = []
    for r in rows("ads/ads_table.csv"):
        if cut and r["cut"] != cut:
            continue
        if mol and r["mol"] != mol:
            continue
        if amorphous_only and r["cell"].startswith("crystal"):
            continue
        r = dict(r); r["Em"] = F(r["Em"]); r["Eads"] = F(r["Eads"]); r["dZ"] = F(r["dZ"])
        r["reg"] = region(r["cell"])
        out.append(r)
    return out


def save(fig, name, note=""):
    import matplotlib.pyplot as plt
    os.makedirs("figs", exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"figs/{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    made.append(f"{name}  {note}")


def fig1():
    """The three patch types side by side: what a "patch" is.

    Top row: side view (x-z), every atom, coloured by element, far atoms faded, so the slab, its free
    top surface and the vacuum above it are visible. Bottom row: top view (x-y) of the exposed set
    only - atoms within DEPTH_EXPOSED of the highest - which is the layer the composition model reads.
    Adsorption sites from ads/sites_<cell>.csv are ringed in both.

    Geometry comes from qe_inputs/<cell>_wf.in, the same relaxed coordinates the DFT used, so nothing
    is recomputed here. Orthographic projection with depth-faded circles: not a ray-traced render, but
    honest and legible at print size.
    """
    import matplotlib.pyplot as plt
    CELLS = os.environ.get("FIG1_CELLS", "crystal_ideal iface_N_D28 matrix_xp").split()
    RAD = {"Mg": 1.60, "Zn": 1.34, "Ca": 1.97}
    DEPTH_EXPOSED = 1.5

    def geom(cell):
        fn = f"qe_inputs/{cell}_wf.in"
        if not os.path.exists(fn):
            return None
        s = open(fn).read()
        cm = re.search(r"CELL_PARAMETERS[^\n]*\n((?:\s*[-\d.eE+]+\s+[-\d.eE+]+\s+[-\d.eE+]+\s*\n){3})", s)
        cell_v = np.array([[float(x) for x in l.split()] for l in cm.group(1).strip().splitlines()]) if cm else None
        blk = s.split("ATOMIC_POSITIONS angstrom", 1)[1].split("K_POINTS", 1)[0]
        el, xyz = [], []
        for l in blk.strip().splitlines():
            q = l.split()
            if len(q) >= 4 and q[0] in RAD:
                el.append(q[0]); xyz.append([float(q[1]), float(q[2]), float(q[3])])
        return np.array(el), np.array(xyz), cell_v

    got = [(c,) + tuple(g) for c in CELLS if (g := geom(c))]
    if not got:
        skipped.append(f"Fig1: no qe_inputs/<cell>_wf.in for {', '.join(CELLS)}"); return
    fig, axes = plt.subplots(2, len(got), figsize=(2.5 * len(got) + 0.5, 5.0),
                             gridspec_kw={"height_ratios": [1.35, 1]})
    axes = np.atleast_2d(axes).reshape(2, len(got))
    for j, (cell, el, P, cv) in enumerate(got):
        sites = set()
        sf = f"ads/sites_{cell}.csv"
        if os.path.exists(sf):
            for r in rows(sf):
                try:
                    sites.add(int(r["atom"]) - 1
                              )
                except (TypeError, ValueError):
                    pass
        # ---- side view, painter's algorithm on y
        ax = axes[0, j]
        order = np.argsort(P[:, 1])
        ymin, ymax = P[:, 1].min(), P[:, 1].max()
        for i in order:
            dep = (P[i, 1] - ymin) / max(ymax - ymin, 1e-9)
            ax.plot(P[i, 0], P[i, 2], "o", ms=3.0 * RAD[el[i]], mfc=C[el[i]], mec="#333333",
                    mew=0.25, alpha=0.30 + 0.70 * dep, zorder=1 + dep)
        for i in sites:
            if i < len(P):
                ax.plot(P[i, 0], P[i, 2], "o", ms=3.0 * RAD[el[i]] + 5, mfc="none",
                        mec="#b2182b", mew=1.3, zorder=50)
        ax.set_title(cell.replace("_", " "), fontsize=8)
        ax.set_aspect("equal"); ax.set_xticks([])
        if j == 0:
            ax.set_ylabel("z (\u00c5)")
        else:
            ax.set_yticklabels([])
        ax.set_ylim(P[:, 2].min() - 1.5, P[:, 2].max() + 5.0)
        ax.text(0.03, 0.97, "vacuum", transform=ax.transAxes, fontsize=6, va="top", color="#888888")
        # ---- top view of the exposed set
        ax = axes[1, j]
        m = (P[:, 2].max() - P[:, 2] < DEPTH_EXPOSED) & (P[:, 2] - P[:, 2].min() >= 2.0)
        for i in np.argsort(P[:, 2]):
            if not m[i]:
                continue
            ax.plot(P[i, 0], P[i, 1], "o", ms=4.2 * RAD[el[i]], mfc=C[el[i]], mec="#333333", mew=0.3)
        for i in sites:
            if i < len(P) and m[i]:
                ax.plot(P[i, 0], P[i, 1], "o", ms=4.2 * RAD[el[i]] + 5, mfc="none",
                        mec="#b2182b", mew=1.3, zorder=50)
        n = int(m.sum())
        frac = " ".join(f"{e} {100*(el[m]==e).mean():.0f}%" for e in EL if (el[m] == e).any())
        ax.set_title(f"outermost {DEPTH_EXPOSED} \u00c5: {n} atoms\n{frac}", fontsize=6.5)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        if cv is not None:
            ax.set_xlim(-1, cv[0][0] + 1); ax.set_ylim(-1, cv[1][1] + 1)
    h = [plt.Line2D([], [], ls="", marker="o", mfc=C[e], mec="#333333", ms=5, label=e) for e in EL]
    h.append(plt.Line2D([], [], ls="", marker="o", mfc="none", mec="#b2182b", mew=1.3, ms=7,
                        label="adsorption site"))
    fig.legend(handles=h, fontsize=7, frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.02))
    save(fig, "fig1_patch_types", f"{len(got)} cells: " + ", ".join(c for c, *_ in got))


def fig2():
    """Mg p-band centre against Zn/Ca neighbour count, AT FIXED DEPTH CLASS, one panel per class.

    The depth stratification is not optional. Depth and neighbour count are correlated - surface atoms
    are under-coordinated - so pooling across depth folds the depth effect into the neighbour slope,
    dilutes it, and in interphase cells reverses its sign. Classes and exclusions follow
    06_analyze_pdos.py exactly: surface < 3 A, sub-surface 3-6 A, interior >= 6 A, with `edge`
    (periodic-seam) and `bottom` (artificial lower free surface) atoms dropped.
    The note printed by `make` lists per-cell within-class slopes to check against the manuscript's
    -0.125 +- 0.010 eV/neighbour for the matrix cells and -0.093 +- 0.040 for the interphase.
    """
    import matplotlib.pyplot as plt
    # e_p_occ, NOT e_p_win. 06_analyze_pdos.py defines e_p_occ as the p-band centre of the OCCUPIED
    # part (first moment from the band bottom to E_F) - the manuscript's "p-band centre", which sits
    # ~2.6-3.0 eV below E_F. e_p_win is the separate "p-window centre" over [E_F-3, E_F+1.5], runs
    # -0.8 to +0.45 eV, and its own docstring warns it is truncated when the SCF lacks bands.
    PB = os.environ.get("PBAND", "e_p_occ")
    SURF, SUB = 3.0, 6.0
    # NOTE: bitwise &, not a chained comparison. "SURF <= d < SUB" expands to
    # "(SURF <= d) and (d < SUB)", and `and` on a numpy array raises ValueError.
    CLS = [("surface (<3 \u00c5)", lambda d: d < SURF),
           ("sub-surface (3\u20136 \u00c5)", lambda d: (d >= SURF) & (d < SUB)),
           ("interior (\u22656 \u00c5)", lambda d: d >= SUB)]
    cells = {}
    for fn in sorted(glob.glob("analysis_*.csv")):
        cell = os.path.basename(fn)[9:-4]
        if region(cell) not in ("interphase", "glass matrix"):
            continue
        rr = rows(fn)
        if not rr or not {"n_ZnCa", PB, "depth"} <= set(rr[0]):
            continue
        pts = []
        for r in rr:
            if r.get("el") != "Mg":
                continue
            if F(r.get("edge", 0)) == 1 or F(r.get("bottom", 0)) == 1:
                continue            # seam atoms and the artificial bottom surface
            v = (F(r["n_ZnCa"]), F(r[PB]), F(r["depth"]))
            if np.isfinite(v).all():
                pts.append(v)
        if len(pts) >= 10:
            cells[cell] = np.array(pts)
    if not cells:
        skipped.append(f"Fig2: no analysis_*.csv has Mg n_ZnCa, {PB} and depth"); return
    RC = {"interphase": "#7570b3", "glass matrix": "#1b9e77"}

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.8), sharey=True)
    for ax, (lab, f) in zip(axes, CLS):
        for r in ("interphase", "glass matrix"):
            P = [p[f(p[:, 2])] for c, p in cells.items() if region(c) == r]
            P = np.vstack([q for q in P if len(q)]) if any(len(q) for q in P) else np.empty((0, 3))
            if len(P) < 8:
                continue
            ax.plot(P[:, 0] + (0.07 if r == "interphase" else -0.07), P[:, 1], ".", ms=2,
                    alpha=0.18, color=RC[r], zorder=1)
            xs, mu, se = [], [], []
            for v in sorted(set(P[:, 0])):
                y = P[P[:, 0] == v, 1]
                if len(y) >= 5:
                    xs.append(v); mu.append(y.mean()); se.append(y.std(ddof=1) / np.sqrt(len(y)))
            if xs:
                ax.errorbar(xs, mu, yerr=se, fmt="o", ms=3.6, capsize=2, lw=1.1,
                            color=RC[r], zorder=3)
            if len(set(P[:, 0])) >= 3:
                k, b = np.polyfit(P[:, 0], P[:, 1], 1)
                xr = np.array([P[:, 0].min(), P[:, 0].max()])
                ax.plot(xr, k * xr + b, "-", lw=1.5, color=RC[r], zorder=4,
                        label=f"{r}  {k:+.3f}")
        ax.set_title(lab, fontsize=8)
        ax.set_xlabel("Zn/Ca neighbours")
        if ax.get_legend_handles_labels()[1]:
            ax.legend(fontsize=6, frameon=False, loc="lower left")
    axes[0].set_ylabel("Mg p-band centre (eV)")

    # per-cell slope with depth-class fixed effects: centre within each class, then one regression
    note, byreg = f"{len(cells)} cells, {sum(len(p) for p in cells.values())} Mg atoms "\
                  f"(edge and bottom atoms excluded), {PB} at fixed depth class", {}
    for c, p in sorted(cells.items()):
        y = p[:, 1].copy()
        for _, f in CLS:
            m = f(p[:, 2])
            if m.sum() >= 3:
                y[m] -= y[m].mean()
        if len(set(p[:, 0])) >= 3:
            byreg.setdefault(region(c), []).append((c, np.polyfit(p[:, 0], y, 1)[0]))
    for r, v in sorted(byreg.items()):
        k = np.array([s for _, s in v])
        note += (f"\n      {r}: {k.mean():+.3f} +- {k.std(ddof=1) / np.sqrt(len(k)):.3f} eV/neighbour  ("
                 + ", ".join(f"{s:+.3f}" for _, s in v) + ")")
    save(fig, "fig2_pband_vs_neighbours", note)


def fig3():
    """Dipole-corrected work function per cell, grouped by region."""
    import matplotlib.pyplot as plt
    vals, src = {}, ""
    if os.path.exists("figs/wf.csv"):
        for r in rows("figs/wf.csv"):
            vals[r["cell"]] = F(r["W"])
        src = "figs/wf.csv"
    else:
        # EVERY work-function log, oldest first, so a later value supersedes an earlier one.
        # Any single log is usually a partial run: the manuscript's claim spans ten cells and
        # no one log holds them all.
        logs = sorted({f for p in ("log.wf_*.out", "log.workfunction_*.out", "log.analyze_*.out",
                                   "log.ads_analyze_*.out", "log.*wf*.out") for f in glob.glob(p)},
                      key=os.path.getmtime)
        for lg in logs:
            cur = None
            for l in open(lg, errors="ignore"):
                m = re.match(r"\s*(\S+): E_F = ", l)
                if m:
                    cur = m.group(1)
                m2 = re.search(r"potential read as Ry.*->\s*W =\s*([-\d.]+) eV", l)
                if m2 and cur:
                    vals[cur.replace("_wf", "")] = float(m2.group(1))
        src = f"{len(logs)} log(s)"
        if vals:
            os.makedirs("figs", exist_ok=True)
            with open("figs/wf.csv", "w", newline="") as f:
                w = csv.writer(f); w.writerow(["cell", "W"])
                for c in sorted(vals):
                    w.writerow([c, f"{vals[c]:.4f}"])
    # name the finished slabs that have no work function anywhere, so the figure cannot quietly
    # under-represent the claim it is supporting
    miss = {os.path.basename(f)[:-7] for f in glob.glob("qe_inputs/*_wf.out")
            if "JOB DONE" in open(f, errors="ignore").read()} - set(vals)
    if miss:
        skipped.append(f"Fig3: {len(miss)} finished slab(s) have no work function in any log: "
                       + ", ".join(sorted(miss))
                       + "\n      -> rerun 10_workfunction.py for those (needs wf/<name>.cube), or add rows to figs/wf.csv")
    if not vals:
        skipped.append("Fig3: no figs/wf.csv and no work-function log found. Produce one with\n"
                       "        python3 10_workfunction.py <name> ... > log.wf_manual.out\n"
                       "      (needs wf/<name>.cube from pp.x plot_num=11), then rerun."); return
    order = sorted(vals, key=lambda c: (("crystal", "interphase", "glass matrix", "series glass").index(region(c)), c))
    fig, ax = plt.subplots(figsize=(max(3.4, 0.42 * len(order) + 1.4), 2.8))
    for i, c in enumerate(order):
        ax.bar(i, vals[c], color={"interphase": "#7570b3", "glass matrix": "#1b9e77",
                                  "crystal": "#999999", "series glass": "#1b9e77"}[region(c)], width=0.68)
    cry = [vals[c] for c in vals if region(c) == "crystal"]
    if cry:
        ax.axhline(np.mean(cry), ls="--", lw=1, color="k")
        ax.text(len(order) - 0.4, np.mean(cry), f" Mg(0001) {np.mean(cry):.2f} eV", fontsize=7, va="bottom", ha="right")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([c.replace("_", " ") for c in order], rotation=60, ha="right", fontsize=6.5)
    ax.set_ylabel("work function (eV)")
    lo = min(vals.values()); hi = max(vals.values())
    ax.set_ylim(max(0, lo - 0.35), hi + 0.35)
    byreg = {}
    for c, v in vals.items():
        byreg.setdefault(region(c), []).append(v)
    note = " | ".join(f"{k} mean {np.mean(v):.2f} eV (n={len(v)})" for k, v in sorted(byreg.items()))
    am = [v for c, v in vals.items() if not c.startswith("crystal")]
    if am:
        note += f" | amorphous span {min(am):.2f}-{max(am):.2f} eV"
    save(fig, "fig3_workfunction", f"{len(order)} cells from {src}\n      {note}")


def fig4():
    """Site-matched H2O binding, all amorphous sites, grouped by element, region distinguished."""
    import matplotlib.pyplot as plt
    d = [r for r in ads("H2O") if np.isfinite(r["Em"])]
    if len(d) < 6:
        skipped.append(f"Fig4: only {len(d)} amorphous H2O sites with Em"); return
    rng = np.random.default_rng(1)
    lab = {}
    fig, ax = plt.subplots(figsize=(3.5, 3.2))
    for i, e in enumerate(EL):
        v = [r for r in d if r["el"] == e]
        if not v:
            continue
        for r in v:
            phys = F(r.get("cn")) == 0          # no metal neighbour within 2.5 A: physisorbed
            ax.plot(i + rng.uniform(-0.17, 0.17), r["Em"], MK[r["reg"]], ms=4.2,
                    mfc=(C[e] if phys else "none"), mec=C[e], mew=1.0, alpha=0.9)
        m = np.mean([r["Em"] for r in v]); s = np.std([r["Em"] for r in v], ddof=1) if len(v) > 1 else 0
        ax.plot([i - 0.32, i + 0.32], [m, m], "-", lw=2.0, color=C[e], solid_capstyle="butt")
        ax.errorbar(i, m, yerr=s, color=C[e], capsize=3, lw=1.2)
        nph = sum(1 for r in v if F(r.get("cn")) == 0)
        lab[i] = (e, len(v), nph)
    ax.axhline(0, color="#bbbbbb", lw=0.8, zorder=0)
    import matplotlib.transforms as mtr
    tr = mtr.blended_transform_factory(ax.transData, ax.transAxes)   # all labels on one line, clear of the data
    for i, (e, k, nph) in lab.items():
        ax.text(i, 1.015, f"n={k}" + (f"\n{nph} physisorbed" if nph else ""), transform=tr,
                ha="center", va="bottom", fontsize=6.5, color=C[e], linespacing=1.25)
    ax.set_xticks(range(len(EL))); ax.set_xticklabels(EL)
    ax.set_xlabel("element under the adsorbate")
    ax.set_ylabel(r"$E_\mathrm{ads}$(H$_2$O), site-matched (eV)")
    h = [plt.Line2D([], [], ls="", marker=MK[k], mfc="none", mec="k", label=k)
         for k in ("interphase", "glass matrix") if any(r["reg"] == k for r in d)]
    if any(F(r.get("cn")) == 0 for r in d):
        h.append(plt.Line2D([], [], ls="", marker="o", mfc="k", mec="k", label="physisorbed"))
    ax.legend(handles=h, fontsize=7, frameon=False, loc="lower right")
    save(fig, "fig4_binding_by_element", f"{len(d)} sites")


def fig5():
    """E_ads against clean-slab excess charge dZ, with the within-Mg fit."""
    import matplotlib.pyplot as plt
    d = [r for r in ads("H2O") if np.isfinite(r["Em"]) and np.isfinite(r["dZ"])]
    if len(d) < 6:
        skipped.append(f"Fig5: only {len(d)} H2O sites have both Em and dZ"); return
    fig, ax = plt.subplots(figsize=(3.5, 2.9))
    for e in EL:
        v = [r for r in d if r["el"] == e]
        if v:
            ax.plot([r["dZ"] for r in v], [r["Em"] for r in v], MK["glass matrix"], ms=4.2,
                    mfc="none", mec=C[e], mew=1.0, label=f"{e} (n={len(v)})")
    mg = [r for r in d if r["el"] == "Mg"]
    if len(mg) >= 4:
        x = np.array([r["dZ"] for r in mg]); y = np.array([r["Em"] for r in mg])
        k, b = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 10)
        rho = np.corrcoef(np.argsort(np.argsort(x)), np.argsort(np.argsort(y)))[0, 1]
        ax.plot(xs, k * xs + b, "-", lw=1.2, color=C["Mg"], alpha=0.8)
        ax.text(0.03, 0.04, f"within Mg: slope {k:+.2f} eV/e, Spearman $\\rho$={rho:+.2f}",
                transform=ax.transAxes, fontsize=6.5, color=C["Mg"])
    ax.set_xlabel(r"clean-slab excess charge $\Delta Z$ (e)")
    ax.set_ylabel(r"$E_\mathrm{ads}$(H$_2$O), site-matched (eV)")
    ax.legend(fontsize=7, frameon=False)
    save(fig, "fig5_binding_vs_charge", f"{len(d)} sites")


def fig6():
    """Barrier profiles from the newest barrier-analyze log."""
    import matplotlib.pyplot as plt
    lg = newest("log.barrier_analyze_*.out")
    if not lg:
        skipped.append("Fig6: no log.barrier_analyze_*.out"); return
    scans, cur = {}, None
    for l in open(lg, errors="ignore"):
        m = re.match(r"^(Ca|Mg|Zn)\s+(\S+?):(s\d+)\s", l)
        if m:
            cur = (m.group(1), f"{m.group(2)}:{m.group(3)}"); scans[cur] = []
            continue
        m = re.match(r"\s*lambda=\s*([\d.]+)\s+([+-][\d.]+)\s+(\S.*)$", l)
        if m and cur:
            scans[cur].append((float(m.group(1)), float(m.group(2)), m.group(3).strip()))
    scans = {k: v for k, v in scans.items() if len(v) >= 3}
    if not scans:
        skipped.append(f"Fig6: parsed no profiles from {lg}"); return
    fig, ax = plt.subplots(figsize=(3.6, 2.9))
    for (el, lab), pts in sorted(scans.items()):
        pts.sort()
        ok = [(l, e) for l, e, s in pts if s in ("ok", "reactant", "product")]
        bad = [(l, e) for l, e, s in pts if s not in ("ok", "reactant", "product")]
        if ok:
            ax.plot(*zip(*ok), "-o", ms=3.6, lw=1.3, color=C[el], label=f"{el}  {lab}")
        if bad:
            ax.plot(*zip(*bad), "x", ms=4.5, color=C[el], alpha=0.6)
        for (l1, e1), (l2, e2) in zip(ok, ok[1:]):          # a near-vertical step = branch switch
            if abs(e2 - e1) > 0.08 and (l2 - l1) < 0.05:
                ax.axvspan(l1, l2, color=C[el], alpha=0.12, lw=0, zorder=0)
        interior = [(l, e) for l, e in ok if 0 < l < 1]
        if interior:
            lm, bm = max(interior, key=lambda t: t[1])
            ax.plot([lm], [bm], "o", ms=7, mfc="none", mec=C[el], mew=1.2)
            ax.annotate(f"{max(bm,0):.2f} eV", (lm, bm), textcoords="offset points",
                        xytext=(4, 5), fontsize=6.5, color=C[el])
    ax.axhline(0, color="#bbbbbb", lw=0.8, zorder=0)
    ax.set_xlabel(r"reaction coordinate $\lambda$   (H$_2$O* $\rightarrow$ OH* + H*)")
    ax.set_ylabel(r"$E$ relative to H$_2$O* (eV)")
    ax.legend(fontsize=6.5, frameon=False)
    save(fig, "fig6_barrier_profiles", f"{len(scans)} scans from {os.path.basename(lg)}")


def _model_rows():
    out = []
    for r in rows("surface_model.csv"):
        cell = re.sub(r"_wf$", "", r["cell"])
        out.append({"cell": cell, "n": int(F(r["n"])), "fZn": F(r["fZn"]), "fMg": F(r["fMg"]),
                    "fCa": F(r["fCa"]), "E": F(r["E_A"]), "lo": F(r["lo"]), "hi": F(r["hi"]),
                    "reg": region(cell), "comp": (re.match(r"(c\d+_[a-z]+\d+)_", cell) or [None, None])[1]})
    return out


def fig7():
    """Predicted patch binding against exposed Zn fraction, with measured patch means overlaid."""
    import matplotlib.pyplot as plt
    M = [r for r in _model_rows() if r["n"] >= 7 and np.isfinite(r["E"])]
    if len(M) < 6:
        skipped.append(f"Fig7: only {len(M)} model rows with n>=7"); return
    # Stratified estimator. Sites were chosen round-robin by element, NOT in proportion to each
    # patch's composition, so a plain mean over them is biased toward whichever element is
    # over-sampled. Reweight each element's site mean by that patch's own exposed fraction,
    # renormalised over the elements that actually have a site.
    sites = {}
    for r in ads("H2O"):
        if np.isfinite(r["Em"]):
            sites.setdefault(r["cell"], {}).setdefault(r["el"], []).append(r["Em"])
    frac = {r["cell"]: {"Ca": r["fCa"], "Mg": r["fMg"], "Zn": r["fZn"]} for r in M}
    meas = {}
    for cell, byel in sites.items():
        f = frac.get(cell)
        if not f:
            continue
        w = sum(f.get(e, 0.0) for e in byel)
        if w <= 0:
            continue
        est = sum(f.get(e, 0.0) / w * np.mean(v) for e, v in byel.items())
        var = sum((f.get(e, 0.0) / w) ** 2 * (np.var(v, ddof=1) / len(v))
                  for e, v in byel.items() if len(v) > 1)
        meas[cell] = (est, float(np.sqrt(var)), sum(len(v) for v in byel.values()), len(byel))
    fig, ax = plt.subplots(figsize=(3.6, 2.9))
    # E_pred depends on fCa and fMg too, so it is NOT a function of fZn alone: scatter the patches,
    # never join them with a line. The binary Mg-Zn mixing limit is single-valued and gives the trend.
    for r in M:
        ax.plot([r["fZn"], r["fZn"]], [r["lo"], r["hi"]], "-", lw=0.6, color="#cccccc", zorder=1)
    ax.plot([r["fZn"] for r in M], [r["E"] for r in M], ".", ms=3.4, color="#444444",
            zorder=2, label="patch prediction")
    ael = {}
    for r in ads("H2O"):
        if np.isfinite(r["Em"]):
            ael.setdefault(r["el"], []).append(r["Em"])
    if {"Mg", "Zn"} <= set(ael):
        eM, eZ = np.mean(ael["Mg"]), np.mean(ael["Zn"])
        xs = np.linspace(0, 1, 50)
        ax.plot(xs, (1 - xs) * eM + xs * eZ, "--", lw=1.0, color="#444444",
                label="Mg\u2013Zn mixing limit")
    # standard error of each patch mean, not the within-patch sd: the sd is dominated by the
    # element effect the model is built on, so it tests nothing
    res = []
    for r in M:
        if r["cell"] in meas:
            est, se, ns, nel = meas[r["cell"]]
            ax.errorbar(r["fZn"], est, yerr=se, fmt=MK[r["reg"]], ms=5, mfc="none",
                        mec="#b2182b", ecolor="#b2182b", mew=1.1, capsize=2.5, zorder=5)
            res.append(est - r["E"])
    ax.plot([], [], MK["glass matrix"], mfc="none", mec="#b2182b",
            label="measured, composition-weighted")
    if res:
        ax.text(0.98, 0.04, f"RMSE {np.sqrt(np.mean(np.square(res))):.3f} eV over {len(res)} patches",
                transform=ax.transAxes, ha="right", fontsize=6.5, color="#b2182b")
    ax.set_xlabel("exposed Zn fraction")
    ax.set_ylabel(r"$E_\mathrm{ads}$(H$_2$O) (eV)")
    ax.legend(fontsize=6.5, frameon=False, loc="upper left")
    save(fig, "fig7_model_vs_exposed_zn", f"{len(M)} patches, {len(meas)} measured"
         + (f", stratified-estimator RMSE {np.sqrt(np.mean(np.square(res))):.3f} eV" if res else ""))


def fig8():
    """Composition series: predicted binding against exposed Zn for every composition."""
    import matplotlib.pyplot as plt
    M = [r for r in _model_rows() if r["comp"] and r["n"] >= 7 and np.isfinite(r["E"])]
    S = rows("series_model.csv")
    if not M or not S:
        skipped.append("Fig8: needs series rows in surface_model.csv and series_model.csv"); return
    comps = sorted({r["comp"] for r in M})
    cmap = {c: plt.get_cmap("viridis")(i / max(1, len(comps) - 1)) for i, c in enumerate(comps)}
    fig, ax = plt.subplots(figsize=(3.8, 3.0))
    for c in comps:
        v = [r for r in M if r["comp"] == c]
        ax.plot([r["fZn"] for r in v], [r["E"] for r in v], "o", ms=3.4, mfc="none",
                mec=cmap[c], mew=0.9, alpha=0.85)
    for r in S:
        c = r["composition"]
        ax.errorbar(F(r["fZn"]), F(r["E_pred"]), yerr=F(r["se"]), fmt="D", ms=5.5,
                    color=cmap.get(c, "k"), capsize=3, mec="k", mew=0.6, zorder=6,
                    label=f"{c}  ({r['patches']} patches)")
    ax.set_xlabel("exposed Zn fraction")
    ax.set_ylabel(r"predicted $E_\mathrm{ads}$(H$_2$O) (eV)")
    ax.legend(fontsize=6, frameon=False, loc="best")
    ax.text(0.98, 0.02, "diamonds: composition means\ncircles: individual patches",
            transform=ax.transAxes, fontsize=6, ha="right", va="bottom", color="#444444")
    save(fig, "fig8_composition_series", f"{len(comps)} compositions, {len(M)} patches")


def figS1():
    """MEAM vs DFT surface relaxation: the per-element displacements and the Zn-Mg splitting."""
    import matplotlib.pyplot as plt
    meam, dft, rel = {}, {}, []
    lg = newest("log.relaxprofile_*.out")
    if lg:
        sec = ""
        for l in open(lg, errors="ignore"):
            if "HEIGHT above" in l:
                sec = "h"; continue
            m = re.match(r"\s*(Ca|Mg|Zn)\s+([-\d.]+)\s+([-\d.]+)\s+([+-][\d.]+)\s+([\d.]+)", l)
            if m and sec == "h":
                meam[m.group(1)] = (float(m.group(4)), float(m.group(5)))
    lg2 = newest("log.dftheight_*.out")
    if lg2:
        sec = ""
        for l in open(lg2, errors="ignore"):
            if "mean DFT displacement" in l:
                sec = "d"; continue
            m = re.match(r"\s*(Ca|Mg|Zn)\s+([+-][\d.]+)\s+\+-\s+([\d.]+)", l)
            if m and sec == "d":
                dft[m.group(1)] = (float(m.group(2)), float(m.group(3)))
            m2 = re.match(r"\S+\s+\d+\s+\d+\s+\|.*\s([+-][\d.]+)\s*$", l)
            if m2:
                rel.append(float(m2.group(1)))
    if not meam and not dft:
        skipped.append("FigS1: needs log.relaxprofile_*.out and log.dftheight_*.out"); return
    fig, axes = plt.subplots(1, 2, figsize=(6.2, 2.7))
    ax = axes[0]
    # The two relaxations run in SEQUENCE: cut -> (MEAM) -> MEAM geometry -> (DFT) -> DFT geometry.
    # Drawing them as side-by-side bars would imply two independent measurements of one quantity.
    # Second bar is therefore CUMULATIVE displacement from the cut, = MEAM + DFT.
    w = 0.36
    for i, e in enumerate(EL):
        if e in meam:
            ax.bar(i - w / 2, meam[e][0], w, yerr=meam[e][1], color=C[e], alpha=0.95, capsize=3)
        if e in meam and e in dft:
            tot = meam[e][0] + dft[e][0]
            ax.bar(i + w / 2, tot, w, yerr=np.hypot(meam[e][1], dft[e][1]),
                   color=C[e], alpha=0.45, capsize=3, hatch="///")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(EL))); ax.set_xticklabels(EL)
    ax.set_ylabel(r"displacement from the cut ($\mathrm{\AA}$, + = outward)")
    ax.set_title("solid: after MEAM     hatched: after MEAM then DFT", fontsize=7.5)
    if {"Mg", "Zn"} <= set(meam) and "Zn" in dft and "Mg" in dft:
        sm = meam["Zn"][0] - meam["Mg"][0]
        st = (meam["Zn"][0] + dft["Zn"][0]) - (meam["Mg"][0] + dft["Mg"][0])
        ax.set_ylim(ax.get_ylim()[0] * 1.18, ax.get_ylim()[1])
        ax.text(0.5, 0.025, f"Zn$-$Mg splitting:  MEAM {sm:+.2f} $\\mathrm{{\\AA}}$   "
                           f"after DFT {st:+.2f} $\\mathrm{{\\AA}}$",
                transform=ax.transAxes, ha="center", fontsize=6.5, color="#333333",
                bbox=dict(fc="white", ec="none", pad=1.5))
    ax = axes[1]
    if rel:
        r = np.array([v for v in rel if abs(v) < 1.0])          # drop the collapsed-atom outlier
        ax.hist(r, bins=14, color="#777777", alpha=0.85)
        ax.axvline(np.median(r), color="#b2182b", lw=1.5)
        ax.text(np.median(r), ax.get_ylim()[1] * 0.95, f" median {np.median(r):+.3f}",
                color="#b2182b", fontsize=7, va="top")
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel(r"DFT $\Delta z$(Zn) $-$ $\Delta z$(Mg) per slab ($\mathrm{\AA}$)")
        ax.set_ylabel(f"slabs (n={len(r)})")
    else:
        ax.axis("off"); ax.text(0.5, 0.5, "per-slab values not parsed", ha="center", fontsize=8)
    save(fig, "figS1_meam_vs_dft_relaxation", f"MEAM {len(meam)} / DFT {len(dft)} elements")


def tables():
    os.makedirs("figs", exist_ok=True)
    n = 0
    M = _model_rows()
    if M:
        with open("figs/table1_patch_inventory.csv", "w", newline="") as f:
            w = csv.writer(f); w.writerow(["patch", "region", "composition", "exposed_atoms", "fCa", "fMg", "fZn", "E_pred", "lo", "hi"])
            for r in sorted(M, key=lambda r: (r["reg"], r["cell"])):
                w.writerow([r["cell"], r["reg"], r["comp"] or "", r["n"],
                            f"{r['fCa']:.4f}", f"{r['fMg']:.4f}", f"{r['fZn']:.4f}",
                            f"{r['E']:.4f}", f"{r['lo']:.4f}", f"{r['hi']:.4f}"])
        n += 1
    for mol in ("H2O", "OH"):
        d = [r for r in ads(mol) if np.isfinite(r["Em"])]
        if not d:
            continue
        with open(f"figs/table2_binding_by_element_{mol}.csv", "w", newline="") as f:
            w = csv.writer(f); w.writerow(["element", "mean_eV", "sd_eV", "n", "n_no_metal_neighbour"])
            for e in EL:
                v = [r for r in d if r["el"] == e]
                if v:
                    w.writerow([e, f"{np.mean([r['Em'] for r in v]):.4f}",
                                f"{np.std([r['Em'] for r in v], ddof=1):.4f}" if len(v) > 1 else "",
                                len(v), sum(1 for r in v if F(r.get('cn')) == 0)])
        n += 1
    d = [r for r in ads(None, amorphous_only=False) if np.isfinite(r["Em"])]
    if d:
        with open("figs/tableSI_all_sites.csv", "w", newline="") as f:
            w = csv.writer(f); w.writerow(["cell", "region", "site", "mol", "element", "dZ", "coordination",
                                           "depth", "zO", "Eads_unmatched", "Eads_matched"])
            for r in sorted(d, key=lambda r: (r["mol"], r["el"], r["cell"], r["site"])):
                w.writerow([r["cell"], r["reg"], r["site"], r["mol"], r["el"], r.get("dZ", ""),
                            r.get("cn", ""), r.get("depth", ""), r.get("zO", ""),
                            f"{r['Eads']:.4f}" if np.isfinite(r["Eads"]) else "",
                            f"{r['Em']:.4f}"])
        n += 1
    made.append(f"{n} table CSVs in figs/")


def make():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"matplotlib unavailable: {type(e).__name__}: {e}\n"
              f"Run with PYTHONNOUSERSITE=1 so ~/.local's numpy 2.x stops shadowing anaconda's.")
        return
    import matplotlib as mpl
    mpl.rcParams.update({"font.size": 8, "font.family": "serif", "axes.linewidth": 0.8,
                         "xtick.direction": "out", "ytick.direction": "out",
                         "axes.spines.top": False, "axes.spines.right": False,
                         "legend.handlelength": 1.4, "figure.dpi": 150})
    for fn in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, figS1, tables):
        try:
            fn()
        except Exception as e:
            skipped.append(f"{fn.__name__}: {type(e).__name__}: {e}")
    print("=" * 78)
    print(f"MADE ({len(made)})")
    for m in made:
        print("  " + m)
    print(f"\nSKIPPED ({len(skipped)})")
    for s in skipped:
        print("  " + s)
    print("=" * 78)
    print("\nfigs/ holds PDF (for the manuscript) and PNG (to eyeball) of each, plus the table CSVs.")


if __name__ == "__main__":
    c = sys.argv[1] if len(sys.argv) > 1 else ""
    if c == "make":
        make()
    elif c == "inspect":
        print("inspect lives in the previous revision; `make` now reports what it could not find.")
        make()
    else:
        print(__doc__)
