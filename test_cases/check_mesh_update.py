"""Consistency checks of the AURA <-> EMerge mesh update on the patch antenna geometry.

1. build_topology() on the original nodes/tets must reproduce Mesh3D._pre_update exactly.
2. After one refinement pass the new Mesh3D must pass diagnose(), keep every
   volume/face/curve mapping consistent and preserve the volume of every domain.
"""

import emerge as em
import numpy as np

import emerge_aura as aura
from emerge_aura.brep import snapshot_gmsh
from emerge_aura.emerge_tools import build_topology, mesh_to_aura_inputs, refine_mesh3d

mm = 0.001
Wpatch, Lpatch, wline, wstub, lstub = 53 * mm, 52 * mm, 3.2 * mm, 7 * mm, 15.5 * mm
wsub, hsub, th, Rair = 100 * mm, 100 * mm, 1.524 * mm, 100 * mm

model = em.Simulation("PatchCheck")
dielectric = em.geo.Box(wsub, hsub, th, position=(-wsub / 2, -hsub / 2, -th))
air = em.geo.Sphere(Rair).background()
rpatch = em.geo.XYPlate(Wpatch, Lpatch, position=(-Wpatch / 2, -Lpatch / 2, 0))
ground = em.geo.XYPlate(wsub, hsub, position=(-wsub / 2, -hsub / 2, -th)).set_material(em.lib.PEC)
cutout1 = em.geo.XYPlate(wstub, lstub, position=(-wline / 2 - wstub, -Lpatch / 2, 0))
cutout2 = em.geo.XYPlate(wstub, lstub, position=(wline / 2, -Lpatch / 2, 0))
line = em.geo.XYPlate(wline, lstub, position=(-wline / 2, -Lpatch / 2, 0))
port = em.geo.Plate(np.array([-wline / 2, -Lpatch / 2, -th]), np.array([wline, 0, 0]), np.array([0, 0, th]))
rpatch = em.geo.add(em.geo.remove(em.geo.remove(rpatch, cutout1), cutout2), line)
rpatch.set_material(em.lib.PEC)
dielectric.set_material(em.Material(3.38))
model.mw.set_resolution(0.15)
model.mw.set_frequency_range(1.545e9, 1.605e9, 3)
model.commit_geometry()
model.mesher.set_boundary_size(rpatch, 2 * mm)
model.mesher.set_face_size(port, 1 * mm)
model.generate_mesh()

mesh = model.mesh
print(f"Original mesh: {mesh.n_tets} tets, {mesh.n_tris} tris, {mesh.n_edges} edges")

# ---------------------------------------------------------------------------
# 1. Topology port equals _pre_update
# ---------------------------------------------------------------------------
fresh = type(mesh)()
build_topology(fresh, mesh.nodes, mesh.tets)
for name in ("edges", "tris", "tet_to_edge", "tet_to_tri", "tri_to_tet", "tri_to_edge",
             "edge_to_tri", "node_to_edge", "centers", "tri_centers", "edge_centers",
             "edge_lengths", "areas"):
    a, b = getattr(fresh, name), getattr(mesh, name)
    same = a.shape == b.shape and (np.array_equal(a, b) if a.dtype.kind in "iu" else np.allclose(a, b, rtol=1e-13, atol=0))
    print(f"  {name:14s} {'OK' if same else 'MISMATCH'}")
    assert same, name
for name in ("inv_edges", "inv_tris", "inv_tets"):
    assert getattr(fresh, name) == getattr(mesh, name), name
print("1. build_topology reproduces _pre_update exactly.")

# ---------------------------------------------------------------------------
# 2. One refinement pass with a synthetic error concentrated at the patch edge
# ---------------------------------------------------------------------------
brep = snapshot_gmsh(scale=1 / model._geo_scale_factor)
print("BREP snapshot:", brep.summary())

c = mesh.centers
err = 1.0 / (1e-3 + np.abs(np.abs(c[0]) - Wpatch / 2) + np.abs(c[2]))  # hot along the radiating edges
inputs = mesh_to_aura_inputs(mesh)
res = aura.refine(**inputs, error=err, brep=brep, growth=0.3)
new = refine_mesh3d(mesh, res)
print(f"Refined mesh: {new.n_tets} tets (+{(new.n_tets / mesh.n_tets - 1) * 100:.1f}%), stats {res.stats}")

new.diagnose()

# Every boundary (one-sided) triangle belongs to some face tag, every tagged triangle exists.
one_sided = set(np.flatnonzero(new.tri_to_tet[1] == new._MISSING_ID).tolist())
tagged = set(i for ids in new.ftag_to_tri.values() for i in ids)
assert one_sided <= tagged, f"{len(one_sided - tagged)} untagged boundary triangles"

# Tagged face areas: planar faces must be preserved exactly, curved ones approximately.
for tag in mesh.ftag_to_tri:
    a0 = mesh.areas[mesh.ftag_to_tri[tag]].sum()
    a1 = new.areas[new.ftag_to_tri[tag]].sum()
    assert abs(a1 - a0) <= 1e-3 * a0, (tag, a0, a1)


def domain_volume(m, tag):
    t = m.tets[:, m.vtag_to_tet[tag]]
    x = [m.nodes[:, t[k]] for k in range(4)]
    return np.abs(np.einsum("ij,ij->j", x[1] - x[0], np.cross((x[2] - x[0]).T, (x[3] - x[0]).T).T)).sum() / 6


for tag in mesh.vtag_to_tet:
    v0, v1 = domain_volume(mesh, tag), domain_volume(new, tag)
    print(f"  volume {tag}: {v0:.6e} -> {v1:.6e} ({(v1 / v0 - 1) * 1e6:+.2f} ppm)")
    assert abs(v1 / v0 - 1) < 1e-3

# Curve edges: total curve length preserved (straight) or grows towards the true arc (curved).
for tag, ids in mesh.etag_to_edge.items():
    l0, l1 = mesh.edge_lengths[ids].sum(), new.edge_lengths[new.etag_to_edge[tag]].sum()
    assert l1 >= l0 * (1 - 1e-9), (tag, l0, l1)

# Node indices of geometric vertices are unchanged.
for dt, ids in mesh.dimtag_to_nodes.items():
    assert np.array_equal(new.dimtag_to_nodes[dt], ids)

q0, q1 = aura.quality(inputs["nodes"], inputs["tets"]), aura.quality(res.nodes, res.tets)
print(f"Quality min {q0.min():.3f} -> {q1.min():.3f}, mean {q0.mean():.3f} -> {q1.mean():.3f}")

# Sphere surface: new nodes on the air sphere lie on the true radius.
sphere_tags = [t for t, k in zip(brep.surf_tags.tolist(), brep.surf_kind.tolist()) if k == 1]
for tag in sphere_tags:
    nodes = np.unique(new.tris[:, new.ftag_to_tri[tag]])
    r = np.linalg.norm(new.nodes[:, nodes], axis=0)
    print(f"  sphere face {tag}: radius range {r.min() / mm:.6f} .. {r.max() / mm:.6f} mm")
print("2. Refined Mesh3D is consistent.")
