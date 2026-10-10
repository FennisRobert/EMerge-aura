# emerge-aura

**AURA — Adaptive Unstructured Refinement Algorithm** for the [EMerge](https://github.com/FennisRobert/EMerge) full-wave FEM solver.

> [!WARNING]
> **This package is vibe coded.**
> Almost all of the code in this repository, including the Rust core, the Python bindings and the EMerge integration, was written by an AI coding assistant (Claude Code). The author directed it and tested the results on a handful of EMerge models, but has **not** reviewed every line by hand.
> It is experimental. Check your results against EMerge's default GMSH-remeshing AMR (`use_aura=False`) before you trust them, and please open an issue if something looks wrong.

## What it does

EMerge's default adaptive mesh refinement (AMR) regenerates the whole mesh with GMSH on every pass, using a size field derived from the error estimate. AURA instead **refines the existing mesh locally**:

0. **Step 0 (once, before the first solve):** narrow strips such as thin PCB traces are detected purely from the geometry tags. The curve edges along each strip are split until they are at most `strip_ratio` (default 4) strip widths long, and the whole mesh gets one quality pass. Without this, every triangle across a narrow strip is a sliver, whatever the triangulation.
1. The tetrahedra with the highest error are bisected through their longest edge, one at a time, until the mesh has `(1 + growth)` times as many tets. Every split is *conforming*: all tets around the split edge are bisected together, so there are never hanging nodes.
2. New nodes on curved faces and curves (spheres, cylinders, circles, B-splines, …) are projected onto a one-time snapshot of the original OCC geometry, so the mesh converges to the real geometry rather than the initial facets.
3. The new tets, plus a few layers of neighbouring old tets where new and old elements meet, are improved by flips and smoothing:
   - 2-3 flips,
   - edge removal for edges with 3 to 7 tets around them,
   - smoothing of interior nodes,
   - in-plane sliding of nodes on planar geometric faces,
   - surface flips (swapping the diagonal between two triangles of the same planar face).

   A change is only accepted if it improves the worst tet it touches. Curve nodes never move, and face nodes never leave their plane.
4. Periodic boundaries are supported. An edge on a periodic face is always split together with all of its periodic images, so paired faces stay exact translated copies. Node numbering stays consistent across each pair, which EMerge needs for consistent edge orientations.
5. Every geometric tag is inherited exactly: volumes, faces (boundary conditions!), curves and points. EMerge's `Mesh3D` is rebuilt from that, **without calling GMSH**.

### Why

Speed. Regenerating the whole mesh every pass is slow on complex geometries with many edges and faces, such as PCBs. The size-field interpolation of the remeshing method also gets very expensive at high tet counts. AURA's cost per pass grows only with the number of tets and doesn't depend on how complicated the geometry is.

The GMSH remeshing AMR builds a fresh, globally optimised mesh every pass, so it remains the reference for mesh quality. AURA aims to stay close to it while being much faster. On smooth geometries with few edges and low tet counts the speed gain is small: about 1.5× for the mesh-update step at 130k tets on the patch antenna test case. That test case also shows AURA keeping the minimum element quality roughly constant over many passes.

## Installation

```bash
pip install emerge-aura
```

Prebuilt wheels are published for Linux (x86_64, aarch64), macOS (Apple Silicon, Intel) and Windows (x64), for Python 3.9 and newer. On other platforms pip builds from source, which needs a Rust toolchain (<https://rustup.rs>). To install the latest development version:

```bash
pip install git+https://github.com/FennisRobert/EMerge-aura.git
```

It will also install as an EMerge extra: `pip install emerge[aura]`.

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
| `aura_options` | `None` | Dict with `max_depth`, `lepp_ratio`, `improve_sweeps`, `improve_layers`, `improve_initial`, `improve_global`, `strip_ratio`, `strip_max_growth`, `project_boundaries` |

For a callable `growth`, `info` contains `n_tets`, `error` (per-tet error indicator), `convergence` (the last max |ΔS|, or `None`) and `auto` (the value `"auto"` would have used). For example:

```python
model.adaptive_mesh_refinement(growth=lambda step, info: 1.0 if step < 3 else 0.3)
```

Periodic models (e.g. `em.HexCell` unit cells) work as well. On `test_cases/test_4_periodic.py` (hexagonal waveguide array cell with a PEC sphere and iris), the AURA AMR result is within 0.15 dB of a uniformly fine reference mesh.

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
src/improve.rs                     Flips (2-3, edge removal 3..7) and constrained smoothing, incl. in-plane sliding
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
python test_cases/check_periodic.py              # periodic face matching / ordering checks (requires EMerge)
```

`check_mesh_update.py` verifies that the GMSH-free topology rebuild reproduces EMerge's own `Mesh3D._pre_update` exactly. It also checks that a refined mesh passes `Mesh3D.diagnose()` and keeps every domain volume, face area and tag intact.

## Releasing

Wheels for all platforms are built by GitHub Actions (`.github/workflows/release.yml`) on every push and pull request. To publish a release to PyPI:

1. Bump `version` in `Cargo.toml`. It is the single source of the package version.
2. Commit, then tag and push: `git tag v0.1.1 && git push origin main v0.1.1`.

The workflow checks that the tag matches `Cargo.toml`, builds and smoke-tests all wheels, and uploads them with the sdist to PyPI using Trusted Publishing.

## Known limitations

- Refinement only: the mesh is never coarsened.
- Quality improvement only covers the refined region plus `improve_layers` rings of neighbours.
- Flat tets in thin layers (e.g. a substrate meshed one element thick) cannot be fixed without refining the whole layer, so they stay as they are.
- Step 0 does not grade the faces next to a refined strip. Very small `strip_ratio` values (around 2) therefore create fan-shaped slivers beside the strip; the default of 4 avoids this.
- The number of tets per pass slightly overshoots the target when the last split has a large shell.

## License

GPL-2.0-or-later, the same as EMerge.
