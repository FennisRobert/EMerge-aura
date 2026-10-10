import emerge as em

"""DEMO 7: Periodic Cells

EMerge supports periodic environments. In constrast to other software, it automates the setup.
In this demonstration we will look at seting up a rectangular waveguide
array in a flat hexagonal periodic Tiling. A section of such an array is
shown below.

+-----------------+-----------------+-----------------+
|  +-----------+  |  +-----------+  |  +-----------+  |
|  |           |  |  |           |  |  |           |  |
|  +-----------+  |  +-----------+  |  +-----------+  |
+--------+--------+--------+--------+--------+--------+
         |  +-----------+  |  +-----------+  |
         |  |           |  |  |           |  |
         |  +-----------+  |  +-----------+  |
+--------+--------+--------+--------+--------+--------+
|  +-----------+  |  +-----------+  |  +-----------+  |
|  |           |  |  |           |  |  |           |  |
|  +-----------+  |  +-----------+  |  +-----------+  |
+--------+--------+--------+--------+--------+--------+

"""
mm = 0.001
a = 106 * mm  # Cell width
b = 30 * mm  # Cell depth
H = 70 * mm  # Airbox Height
wga = 70 * mm  # Waveguide width
wgb = 18 * mm  # Waveguide Height
fl = 25 * mm  # Feed length

# We start again by defining our simulation model
model = em.Simulation("Periodic")
model.check_version("3.0.0")  # Checks version compatibility.

periodic_cell = em.HexCell((-a / 2, b / 2, 0), (-a / 2, -b / 2, 0), (0, -b / 2, 0))

model.set_periodic_cell(periodic_cell)

box = periodic_cell.volume(0, H)
waveguide = em.geo.Box(wga, wgb, fl, (-wga / 2, -wgb / 2, -fl))

# Small PEC features that stay clear of the periodic walls, to give the AMR
# something to refine: a metal sphere above the aperture (curved surface, tests
# BREP projection) and a capacitive iris made of two thin plates in the waveguide.
w_iris = 12 * mm
sphere = em.geo.Sphere(6 * mm, (0, 0, 25 * mm)).set_material(em.lib.PEC)
iris_l = em.geo.XYPlate(w_iris, wgb, position=(-wga / 2, -wgb / 2, -fl / 2)).set_material(em.lib.PEC)
iris_r = em.geo.XYPlate(w_iris, wgb, position=(wga / 2 - w_iris, -wgb / 2, -fl / 2)).set_material(em.lib.PEC)
model.commit_geometry()

model.mw.set_frequency_range(2.8e9, 3.3e9, 5)
model.mw.set_resolution(0.1)

# Then we create our mesh and view the result
model.generate_mesh()

# ---------------------------------------------------------------------------
# Periodic consistency checks
# ---------------------------------------------------------------------------
import numpy as np
from emerge._emerge.mth.pairing import pair_coordinates
from emerge_aura.brep import snapshot_gmsh
from emerge_aura.emerge_tools import periodic_report, refine_mesh_step

PBC_DSMAX = 1e-15  # tolerance used by EMerge's assembler


def emerge_pairing_ok(mesh, periodic):
    """Pair tris and edges exactly like EMerge's assembler and verify every pair."""
    ok = True
    for tags1, tags2, dv in periodic:
        for centers, ids1, ids2 in (
            (mesh.tri_centers, mesh.get_triangles(tags1), mesh.get_triangles(tags2)),
            (mesh.edge_centers, mesh.get_edges(tags1), mesh.get_edges(tags2)),
        ):
            m = pair_coordinates(centers, ids1, ids2, np.asarray(dv), PBC_DSMAX)
            d = np.array([np.linalg.norm(centers[:, i] + dv - centers[:, j]) for i, j in m.items()])
            if len(ids1) != len(ids2) or len(set(m.values())) != len(ids2) or d.max() > PBC_DSMAX:
                print(f"   pairing FAILED: {len(ids1)} vs {len(ids2)}, max dist {d.max():.2e}")
                ok = False
    return ok


periodic = [(bc.face1.tags, bc.face2.tags, np.asarray(bc.dv)) for bc in model.mesher._get_periodic_bcs()]
mesh = model.mesh
print(f"Original mesh: {mesh.n_tets} tets, {len(periodic)} periodic pairs")
for r in periodic_report(mesh, periodic):
    print("  ", r)
assert all(r["order_preserving"] for r in periodic_report(mesh, periodic)), "EMerge mesh not order preserving"
assert emerge_pairing_ok(mesh, periodic)

brep = snapshot_gmsh(scale=1 / model._geo_scale_factor)
for p in range(3):
    # Error concentrated near the periodic walls (|x| large or |y| near b/2) and the sphere.
    c = mesh.centers
    err = 1 / (1e-3 + np.abs(np.abs(c[1]) - b / 2)) + 1 / (1e-3 + np.abs(np.linalg.norm(c - [[0], [0], [25 * mm]], axis=0) - 6 * mm))
    mesh, rep = refine_mesh_step(mesh, err, brep=brep, growth=0.5, periodic=periodic)
    mesh.diagnose()
    reports = periodic_report(mesh, periodic)
    print(f"Pass {p + 1}: {rep['n_tets']} tets, q_min {rep['q_min']:.3f}, splits {rep['n_splits']}, projected {rep['n_projected']}")
    for r in reports:
        print("  ", r)
    assert all(r["unmatched"] == 0 and r["order_preserving"] for r in reports)
    assert emerge_pairing_ok(mesh, periodic)
r = np.linalg.norm(mesh.nodes - np.array([[0], [0], [25 * mm]]), axis=0)
sph = [t for t in mesh.ftag_to_tri if np.allclose(r[np.unique(mesh.tris[:, mesh.ftag_to_tri[t]])], 6 * mm, rtol=1e-6)]
print("PEC sphere faces with all nodes on R=6mm:", sph)
print("Periodic refinement is consistent.")
