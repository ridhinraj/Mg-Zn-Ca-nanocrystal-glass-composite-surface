# Water activation on Mg–Zn–Ca nanocrystal/glass composite surfaces

Scripts, structures and derived data for the study of H<sub>2</sub>O and OH adsorption on
surfaces of a dual-phase (amorphous + nanocrystalline) Mg–Zn–Ca composite, combining
2NN MEAM molecular dynamics for the surface structures with dipole-corrected plane-wave
DFT for the electronic structure and adsorption energetics.

> **Paper:** R.V. Raj, A. Zagabathuni, *Exposed element identity, not the nanocrystal–glass
> interphase, controls water activation on Mg–Zn–Ca composite surfaces* (submitted).
> Citation and DOI will be added on acceptance.

## What this repository supports

The central claim is that water binding on these surfaces is set by the **element directly
beneath the adsorbing molecule**, not by whether the site lies in the nanocrystal, the
nanocrystal/glass interphase or the glass matrix. Everything needed to reproduce that result —
and the composition-series model built on it — is here, except the raw DFT output files
(see *Data not included*).

## Layout

```
scripts/     the full pipeline, numbered in execution order
slurm/       prep_*.sh, run_*.sh, submit_*.sh job scripts (SLURM)
structures/  relaxed surface patches: slab_*.xyz and patches/*.data
qe_inputs/   every Quantum ESPRESSO input file used
data/        per-atom descriptors, the per-site adsorption table, model outputs
figures/     the manuscript figures (PDF + PNG) and the table CSVs
```

## Pipeline

Run in numerical order. Scripts 15 and 17 do not exist; the numbering has gaps where
exploratory steps were dropped.

| Step | Script | What it does |
|---|---|---|
| 01–02 | `01_*.in`, `02_relax_surface.in` | cut the equatorial slab from the parent composite, insert vacuum, relax (LAMMPS) |
| 03–04 | `03_region_select_numpy.py`, `04_xyz_to_qe.py` | identify crystal / interphase / glass regions; convert to QE format |
| 05 | `05_make_qe_inputs.py` | SCF and projwfc inputs |
| 06 | `06_analyze_pdos.py` | per-atom Löwdin charges, p-band centres, neighbour counts → `analysis_<cell>.csv` |
| 07–08 | `07_build_periodic_patch.py`, `08_relax_patch.in` | cut laterally periodic patches, relax with MEAM |
| 09–10 | `09_make_wf_inputs.py`, `10_workfunction.py` | dipole-corrected slabs; work function from the vacuum potential |
| 11 | `11_adsorption.py` | site selection, H<sub>2</sub>O/OH placement, **site-matched clean references**, dissociation inputs, analysis → `data/ads_table.csv` |
| 12 | `12_surface_model.py` | the exposed-composition model and its validation |
| 13 | `13_barrier.py` | frozen-hydrogen relaxed scans for the O–H scission barrier |
| 14 | `14_coordination.py` | descriptor tests against the element baseline |
| 16 | `16_checks.py` | k-point, cutoff and PBE+D3 convergence checks |
| 18–20 | `18_segtest.py`, `19_relaxprofile.py`, `20_dftheight.py` | is the outermost-layer Zn excess made by the cut or the relaxation, and does DFT reproduce it |
| 21, 24 | `21_oh.py`, `24_oh_mode.py` | OH on every H<sub>2</sub>O site; element vs binding-geometry |
| 22 | `22_figures.py` | all manuscript figures and tables |
| 23, 25 | `23_pband_region.py`, `25_region_cluster.py` | region comparisons on independent realisations rather than sites |

## Two methodological points worth knowing

**Site-matched clean references.** Adsorption energies are referenced to a clean slab relaxed
with *exactly the same free atoms* as the adsorbate calculation. Using an unmatched reference
shifts amorphous-site energies by 0.19–0.45 eV, which is larger than the effect being measured.
`11_adsorption.py cleanref` builds these and `ads/cref_map.csv` records the mapping.

**Do not re-run `cleanref` on cells that already have computed references.** It names references
`cref_<cell>_<k>` by position in a sorted glob, so adding inputs renumbers them and silently
invalidates finished energies. `21_oh.py` mints new references additively as `cref_<cell>_oh<k>`
for this reason.

## Key data files

| File | Contents |
|---|---|
| `data/ads_table.csv` | one row per adsorption calculation: `cell, site, mol, cut, Eads, Em, el, dZ, n_ZnCa, zO, cn, nh`. `Em` is the site-matched energy used throughout the paper |
| `data/analysis_<cell>.csv` | per atom: position, depth, Löwdin charge, `e_p_occ` (occupied p-band centre), `e_p_win`, Zn/Ca neighbour count |
| `data/surface_model.csv` | per patch: exposed fractions and predicted binding with bootstrap interval |
| `data/series_model.csv` | per composition: pooled exposed fractions and prediction |
| `data/sites_<cell>.csv` | the adsorption sites, with the element and descriptors of each |

`e_p_occ` is the quantity the paper calls the p-band centre (first moment of the occupied
p-PDOS, ≈ −2.6 to −3.0 eV); `e_p_win` is a separate windowed measure and can be truncated when
the SCF lacks bands. Analyses of the p-band use `e_p_occ`.

## Software

- Quantum ESPRESSO 7.5 — PBE, SSSP 1.3.0 efficiency pseudopotentials, 40 Ry / 200 Ry cutoffs,
  Γ-only, Marzari–Vanderbilt smearing, dipole correction on all slab calculations
- LAMMPS (22 Jul 2025) with the 2NN MEAM Mg–Zn–Ca potential of Jang, Seol & Lee,
  *Calphad* **67** (2019) 101674, taken from the NIST Interatomic Potentials Repository,
  record `2019--Jang-H-S-Seol-D-Lee-B-J--Mg-Zn-Ca`
- Python 3 with NumPy and Matplotlib only. No SciPy: all statistics are permutation and
  bootstrap tests written against NumPy

The analysis scripts need NumPy 1.x, because the Matplotlib build they run against is compiled
for it. On a system where a NumPy 2.x install shadows it, run with `PYTHONNOUSERSITE=1`.

## Data not included

Raw Quantum ESPRESSO output files (several GB) and the parent composite structure
(~52 MB, produced by the companion mechanical study) are not in this repository. Every energy
extracted from those outputs is in `data/ads_table.csv`, and the inputs that generated them are
in `qe_inputs/`. The raw outputs are available from the corresponding author on request.

## Reproducing the figures

The scripts were written to run in a flat working directory (`ads/ads_table.csv`, `analysis_*.csv`
and the rest at the top level). This repository groups those files under `data/` for readability,
so build the expected layout once with symlinks:

```bash
bash setup_workdir.sh
cd work
PYTHONNOUSERSITE=1 python3 ../scripts/22_figures.py make
```

Figures are written to `work/figs/`. Everything needed is in the repository, including
`figs/wf.csv` for Figure 3.

### Re-running the DFT

The files in `qe_inputs/` carry the absolute `pseudo_dir` path of the machine they were generated
on. Point them at your own SSSP installation before use:

```bash
sed -i "s|pseudo_dir = '.*'|pseudo_dir = '/path/to/your/pseudo'|" qe_inputs/*.in
```

Adsorption energies are referenced to the clean slabs listed in `data/ads/cref_map.csv`; running an
adsorption input without its matched reference will not reproduce the published numbers.

## Licence

Code (everything in `scripts/` and `slurm/`) is released under the MIT Licence; see `LICENSE`.
Data, structures and figures are released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/): reuse freely with attribution to the
paper above.

## Citation

Please cite the paper above. The interatomic potential should be cited as Jang, Seol & Lee
(2019) and Jang, Kim & Lee, *Calphad* **60** (2018) 200–207.

## Contact

Ridhin V Raj and Aparna Zagabathuni, Department of Materials Science and Engineering,
National Institute of Technology Calicut, Calicut, Kerala 673601, India.
Corresponding author: aparnazagabathuni@nitc.ac.in
