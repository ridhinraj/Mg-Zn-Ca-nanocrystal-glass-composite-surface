"""
Stage 1d (handoff to Stage 2): convert each exported region .xyz file into a
Quantum ESPRESSO pw.x input skeleton, with the bottom ~10 A of atoms tagged
as fixed (matching the same bulk-constraint convention used in the LAMMPS
surface relaxation, step 02).

This produces *skeleton* QE inputs only - ecutwfc/ecutrho, smearing, and the
exact pseudopotential filenames still need to be finalized against your
actual QE pseudopotential directory (SSSP efficiency library recommended for
Mg, Zn, Ca, Sr) in Stage 2. Kept deliberately minimal here so Stage 1 has a
clean, testable output before we tune DFT parameters.

Requires: ase (pip install ase --break-system-packages)
"""

import glob
from ase.io import read
from ase.constraints import FixAtoms

FIXED_THICKNESS = 6.0    # A; patches are ~16 A deep, so fix the bottom ~6 A (revisit in Stage 2)
VACUUM_PAD = 15.0        # extra vacuum to add if the exported xyz lost it

# TODO: point this at your actual QE pseudopotential directory
PSEUDO_DIR = "/path/to/qe_pseudopotentials/SSSP_efficiency/"

# TODO: confirm exact filenames present in your PSEUDO_DIR (names vary by
# SSSP version/vendor - list your directory and match these).
# NOTE: AN_composite_relaxed_300K.data is ternary (Mg/Zn/Ca only, no Sr
# confirmed) - only these three need pseudopotentials for this structure.
# Add Sr back in once you're processing the quaternary Tier-2 compositions.
PSEUDOPOTENTIALS = {
    "Mg": "Mg.pbesol-spn-kjpaw_psl.1.0.0.UPF",
    "Zn": "Zn.pbesol-dn-kjpaw_psl.1.0.0.UPF",
    "Ca": "Ca.pbesol-spn-kjpaw_psl.1.0.0.UPF",
}


def write_qe_relax_input(atoms, prefix, outdir="qe_inputs"):
    import os
    os.makedirs(outdir, exist_ok=True)

    z = atoms.positions[:, 2]
    fixed_mask = z < (z.min() + FIXED_THICKNESS)
    atoms.set_constraint(FixAtoms(mask=fixed_mask))
    print(f"{prefix}: {fixed_mask.sum()} fixed / {len(atoms)} total atoms")

    # Ensure enough vacuum remains along z after cropping in step 03
    cell = atoms.get_cell()
    if cell[2, 2] < (z.max() - z.min() + VACUUM_PAD):
        cell[2, 2] = (z.max() - z.min()) + VACUUM_PAD
        atoms.set_cell(cell)

    input_data = {
        "control": {
            "calculation": "relax",
            "prefix": prefix,
            "pseudo_dir": PSEUDO_DIR,
            "outdir": f"./{prefix}_tmp",
            "tstress": True,
            "tprnfor": True,
            "forc_conv_thr": 1.0e-3,
        },
        "system": {
            "ecutwfc": 50,          # TODO: converge in Stage 2
            "ecutrho": 400,         # TODO: converge in Stage 2 (8x ecutwfc for PAW)
            "occupations": "smearing",
            "smearing": "gaussian",
            "degauss": 0.01,
            "assume_isolated": "2D",   # slab dipole correction, exposed surface
        },
        "electrons": {
            "conv_thr": 1.0e-8,
            "mixing_beta": 0.3,
        },
        "ions": {
            "ion_dynamics": "bfgs",
        },
    }

    from ase.io.espresso import write_espresso_in
    fname = f"{outdir}/{prefix}.relax.in"
    write_espresso_in(
        fname,
        atoms,
        input_data=input_data,
        pseudopotentials=PSEUDOPOTENTIALS,
        kpts=(1, 1, 1),   # Gamma-only: defensible for disordered, non-periodic slab
    )
    print(f"  -> wrote {fname}")


if __name__ == "__main__":
    xyz_files = sorted(glob.glob("slab_*_*.xyz"))
    if not xyz_files:
        print("No exported region .xyz files found. Run 03_phase_id_and_region_select.py first.")
    for f in xyz_files:
        prefix = f.replace(".xyz", "")
        atoms = read(f)
        write_qe_relax_input(atoms, prefix)

    print("\nAll available region slabs converted to QE relax inputs under ./qe_inputs/")
    print("Next: verify PSEUDO_DIR and PSEUDOPOTENTIALS filenames, then move to Stage 2")
    print("(ecutwfc/ecutrho convergence, k-point check, PDOS + Bader charge extraction).")
