# emerge-aura

**AURA — Adaptive Unstructured Refinement Algorithm** for the [EMerge](https://github.com/FennisRobert/EMerge) full-wave FEM solver.

> [!WARNING]
> **This package is vibe coded.**
> Almost all of the code in this repository, including the Rust core, the Python bindings and the EMerge integration, was written by an AI coding assistant (Claude Code). The author directed it and tested the results on a handful of EMerge models, but has **not** reviewed every line by hand.
> It is experimental. Check your results against EMerge's default GMSH-remeshing AMR (`use_aura=False`) before you trust them, and please open an issue if something looks wrong.

## What it does

EMerge's default adaptive mesh refinement (AMR) regenerates the whole mesh with GMSH on every pass, using a size field derived from the error estimate. AURA instead **refines the existing mesh locally**:

1. The tetrahedra with the highest error are bisected through their longest edge, one at a time, until the mesh has `(1 + growth)` times as many tets. Every split is *conforming*: all tets around the split edge are bisected together, so there are never hanging nodes.
2. New nodes on curved faces and curves (spheres, cylinders, circles, B-splines, …) are projected onto a one-time snapshot of the original OCC geometry, so the mesh converges to the real geometry rather than the initial facets.
3. The new tets are improved by 2-3, 3-2 and 4-4 flips and by smoothing interior nodes. A change is only accepted if it improves the worst tet it touches.
4. Every geometric tag is inherited exactly: volumes, faces (boundary conditions!), curves and points. EMerge's `Mesh3D` is rebuilt from that, **without calling GMSH**.

On the patch antenna test case it keeps the minimum element quality roughly constant over many passes. At the same tet count it gives more accurate S-parameters than the GMSH remeshing AMR.

## Installation

You need a Rust toolchain (<https://rustup.rs>), at least until prebuilt wheels are published:

```bash
pip install git+https://github.com/FennisRobert/EMerge-aura.git
# or, from a local clone
pip install ./EMerge-aura
```

Once it is on PyPI it will also install as an EMerge extra: `pip install emerge[aura]`.

`emerge-aura` only depends on `numpy`. It does not import EMerge.

## Usage with EMerge

Nothing changes in your script. When `emerge_aura` is installed, `Simulation.adaptive_mesh_refinement` uses it automatically:

```python
model.adaptive_mesh_refinement(
    max_steps=10,
    frequency=1.59e9,
    growth=0.3,          # +30 % tets per pass
    show_mesh=True,      # show the mesh after every pass
)
```

AURA-specific arguments:

| Argument | Default | Meaning |
|---|---|---|
| `growth` | `0.3` | Tet growth per pass: a float, a per-pass list such as `[1.0, 0.6, 0.3]`, `"auto"` (predicted from the error estimate using `theta`/`order`), or a callable `f(step, info) -> float` |
| `growth_limits` | `(0.1, 1.0)` | Clamp on the growth of every pass |
| `use_aura` | `True` | Set to `False` to force the original GMSH remeshing AMR |
| `aura_options` | `None` | Dict with `max_depth`, `lepp_ratio`, `improve_sweeps`, `project_boundaries` |

For a callable `growth`, `info` contains `n_tets`, `error` (per-tet error indicator), `convergence` (the last max |ΔS|, or `None`) and `auto` (the value `"auto"` would have used). For example:

```python
model.adaptive_mesh_refinement(growth=lambda step, info: 1.0 if step < 3 else 0.3)
```

Models with periodic boundary conditions are not supported yet. EMerge automatically falls back to GMSH remeshing for those.

## Standalone API

The core works on plain numpy arrays (row-major: `nodes (N, 3)`, `tets (M, 4)`):

```python
import emerge_aura

result = emerge_aura.refine(
    nodes, tets, error,
    tet_region=region,                              # optional, e.g. volume tag per tet
    face_tris=tris, face_tri_tags=tri_tags,         # optional, tagged boundary triangles
    curve_edges=edges, curve_edge_tags=edge_tags,   # optional, tagged curve edges
    growth=0.3,
)
result.nodes, result.tets, result.tet_region, result.face_tris, result.face_tri_tags, result.stats

q = emerge_aura.quality(result.nodes, result.tets)  # mean-ratio quality per tet (1 = regular)
```

The original nodes keep their indices, and new nodes are appended after them.

## Layout

```
src/refine.rs                      Longest-edge bisection with conforming edge splits and tag inheritance
src/improve.rs                     Flips (2-3, 3-2, 4-4) and constrained smart Laplacian smoothing
src/geom.rs                        BREP snapshot primitives and point projection
src/lib.rs                         PyO3 bindings
python/emerge_aura/__init__.py     refine(), quality()
python/emerge_aura/brep.py         One-time geometry snapshot from GMSH/OCC (gmsh imported lazily)
python/emerge_aura/emerge_tools.py Mesh3D rebuild without GMSH, refine_mesh_step(), growth rules
test_cases/                        EMerge models used for testing, plus check_mesh_update.py
```

## Development

```bash
cargo test --no-default-features                 # Rust unit tests (no Python needed)
pip install .                                    # build the extension with maturin
python test_cases/check_mesh_update.py           # Mesh3D consistency checks (requires EMerge)
```

`check_mesh_update.py` verifies that the GMSH-free topology rebuild reproduces EMerge's own `Mesh3D._pre_update` exactly. It also checks that a refined mesh passes `Mesh3D.diagnose()` and keeps every domain volume, face area and tag intact.

## Known limitations

- Refinement only: the mesh is never coarsened.
- No periodic boundary conditions yet.
- Quality improvement only touches tets created in the current pass.
- The number of tets per pass slightly overshoots the target when the last split has a large shell.

## License

GPL-2.0-or-later, the same as EMerge.
