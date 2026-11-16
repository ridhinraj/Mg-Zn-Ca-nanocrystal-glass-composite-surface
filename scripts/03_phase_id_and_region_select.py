"""
Stage 1c (REVISED): extract three DFT sub-slabs from the single equatorial
slab, using radial distance in the (x,y) plane from the sphere's in-plane
center (103.5, 103.5) - confirmed from the composition analysis of
AN_composite_relaxed_300K.data:

    r <  20 A : pure Mg crystal core
    r = 20-30 A : interphase transition shell
    r >  30 A : bulk glass matrix baseline

Picking crop radii safely inside each zone (not right at the boundary,
to avoid ambiguous atoms):

    crystal_core patch : lateral r < 15 A   (well inside the pure-Mg zone)
    interphase patch   : lateral r = 21-27 A (well inside the transition shell)
    glass_core patch   : lateral r = 45-65 A (safely deep in matrix, still
                          comfortably clear of both the sphere and any
                          concern about periodic self-interaction at ~half
                          the box width)

PTM structure typing is still run and reported as a cross-check - if the
region a crop lands in doesn't match its expected structure type (e.g. the
"crystal_core" patch comes back mostly "Other" rather than HCP), that's a
sign the surface reconstruction during relaxation disordered that region
more than expected, and the crop radii or region definitions should be
revisited before trusting the DFT results.
"""

from ovito.io import import_file, export_file
from ovito.modifiers import (
    PolyhedralTemplateMatchingModifier,
    ExpressionSelectionModifier,
    ComputePropertyModifier,
)

# ---------------------------------------------------------------- SETTINGS
INPUT_DUMP = "equatorial_slab_relaxed.dump"
COMPOSITION_LABEL = "Mg60Zn35Ca5"   # ternary, confirmed - no Sr in this structure

CENTER_X = 103.5
CENTER_Y = 103.5
SURFACE_Z = 103.5          # the new free surface (equatorial cut plane)
SURFACE_WINDOW = 25.0      # keep atoms within this depth of the surface

REGION_RADII = {
    "crystal_core": (0.0, 15.0),
    "interphase":   (21.0, 27.0),
    "glass_core":   (45.0, 65.0),
}

CRYSTAL_TYPE_ID = 2  # PTM HCP id - TODO verify against your OVITO version/GUI

# --------------------------------------------------------------- PIPELINE
pipeline = import_file(INPUT_DUMP)
pipeline.modifiers.append(PolyhedralTemplateMatchingModifier())

# radial distance from the sphere's in-plane center, as a derived property
pipeline.modifiers.append(
    ComputePropertyModifier(
        output_property="RadialDist",
        expressions=[f"sqrt((Position.X-{CENTER_X})^2 + (Position.Y-{CENTER_Y})^2)"],
    )
)

data = pipeline.compute()
zhi = data.cell[2, 3] + data.cell[2, 2]
print(f"Slab surface (equatorial cut) at approx z = {SURFACE_Z}")
print(f"Cell z-top (vacuum boundary) at approx z = {zhi:.2f}")

for name, (rlo, rhi) in REGION_RADII.items():
    expr = (
        f"RadialDist > {rlo} && RadialDist < {rhi} "
        f"&& Position.Z > {SURFACE_Z - SURFACE_WINDOW} && Position.Z <= {SURFACE_Z}"
    )
    pipeline.modifiers.append(ExpressionSelectionModifier(expression=expr))
    frame_data = pipeline.compute()
    n_selected = frame_data.attributes["ExpressionSelection.count"]

    # cross-check: what fraction of the selected patch is PTM-crystalline?
    sel = frame_data.particles["Selection"].array.astype(bool)
    struct = frame_data.particles["Structure Type"].array
    if n_selected > 0:
        crystal_frac = (struct[sel] == CRYSTAL_TYPE_ID).sum() / n_selected
    else:
        crystal_frac = float("nan")

    print(f"\n{name}: r=({rlo}-{rhi} A), {n_selected} atoms selected, "
          f"{100*crystal_frac:.1f}% PTM-crystalline")
    if name == "crystal_core" and crystal_frac < 0.5:
        print("  WARNING: expected mostly crystalline here - check relaxation "
              "or narrow the radius further inward.")
    if name == "glass_core" and crystal_frac > 0.1:
        print("  WARNING: unexpected crystalline content in matrix patch - "
              "check radius is far enough from the sphere.")

    export_file(
        pipeline,
        f"slab_{COMPOSITION_LABEL}_{name}.xyz",
        "xyz",
        columns=["Particle Type", "Position.X", "Position.Y", "Position.Z"],
        only_selected=True,
    )
    pipeline.modifiers.pop()

print("\nDone. Inspect atom counts and crystal_frac above before proceeding")
print("to 04_xyz_to_qe.py. Target ~150-250 atoms per patch; widen/narrow the")
print("radial windows in REGION_RADII if a patch is far off target.")
