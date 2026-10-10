# AURA
Adaptive Unstructured Refinement Algorithm.

This program will function as an adaptive mesh algorithm for the EMerge full wave solver. The current source code is in:

```
/Users/robertfennis/EMerge/EMerge
```

The current version of EMerge uses a complete re-mesh with a desired field size on specific tets. The goal of AURA will be to implement an actual remeshing algorithm on a specific set of tets that moves towards a goal growth rate (+30% but configurable).

# Architecture
The code will be a Rust module, lightweight with as few external dependencies as possible. The rust module should be able to function by taking two datasets from EMerge:

```
 - Mesh object (Mesh3D)
 - Tetrahedra error indicator
 - Goal growth rate in percentage (0.3 = 30%)
```
Aura will then return data that can be used to update the Mesh3D object from EMerge with all derived fields etc.

## Auxilliary mesh data
In how EMerge is designed, all mesh information related to the original geometry is shoved into the Mesh3D object so that GMSH (which originally generated the mesh) can be ditched. Reason being that saving and loading meshes from GMSH does not properly preserve geometric dimension tag (dim, tag) identities. So which boundaries (dim,tag) belong to which triangles etc.

For that reason **all information should be contained in the Mesh3D object**. Likely all required information to do a proper topology preserving refinement is already there. If not it should be extracted in the mesh information update routine: `mesh._pre_update()`. 

Likewise, the updated mesh should be able to also update all the original fields such as the dicts that map geometric entities to mesh entities etc.

Periodic meshes possibly make this more difficult, we may skip that for later or if crucial bake that in from the start.

# Software dependencies
AURA may not import EMerge as a dependency, only as a development dependency during testing. Aura will eventually be uploaded to PyPI or Github only and installed as an optional dependency.

Aura should have a Python API for the Rust module that does the refinement thus.

Mesh3D -> Data -> Aura (Python) -> Aura (Rust) -> AMR code -> Aura (Python) -> Data to update Mesh3D.

There may be an additional function designed that will be implemented into EMerge (as I'm the author) that can take the information from AURA to update all the mesh mappings.

# How I think it may work
Take this as coarse inspiration, you probably know more about this than I do.

As I understand it, the way to correctly build meshes is to first mesh all edges, then all boundaries and then the volume so logically how AMR can work without creating low-quality tets is to: 
 1. Identify bad tet elements
 2. Grow the regions by including N-adjacency layers
 3. Identify which mesh elements (edges and boundaries) are locked on geometric edges through information from the Mesh3D object.
 4. Reconstruct the edge, then boundary and then finally volume mesh by appropriate point insertion.
 5. Compute update information for each region and pass this out of the API for EMerge to consume and update the Mesh3D object.

 # How to develop this
 There is a patch_test.py which you can use to run the AMR algorithm on. Best is probably to monkey-patch certain classes of the original EMerge source code locally so that I can later implement those updates into the Actual EMerge codebase.

** YOU ARE NOT TO TOUCH ANY EMERGE CODE! ONLY READ IT **

Try as much as possible to just monkey patch `Simulation3D` and `Mesh3D` in some way to update the current adaptive mesh refinement algorithm to use the new code and then the mesh update routines required and use patch_test.py to run the tests on. Later on we can incorporate more models.

This folder should contain all the files such that someone with a rust compiler can do `pip install /Users/robertfennis/EMerge/aura` and have it work plug in and ready. We will later install github actions to publish on all operating systems.

# Research Papers
If you need it, I can look for papers on AMR algorithms or other meshing algorithms in order to figure out how to improve local mesh quality without degrading the element quality.