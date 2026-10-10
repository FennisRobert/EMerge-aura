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
model.mw.set_resolution(0.3)

# Then we create our mesh and view the result
model.generate_mesh()
model.view()

wgbc = model.mw.bc.RectangularWaveguide(waveguide.bottom, 1)

# And then the absorbing boundary at the top
abc = model.mw.bc.AbsorbingBoundary(box.back)
periodic_cell.set_scanangle(30, 45)
model.view(bc=True)

model.adaptive_mesh_refinement(minimum_steps=6, show_mesh=False)
data = model.mw.run_sweep()

model.display.add_object(waveguide)
model.display.add_object(box)
model.display.animate().cbar("Ey").add_field(
    data.field[0].cutplane(1 * mm, y=0).scalar("Ey", "complex"), symmetrize=True
)
model.display.animate().cbar("Ey").add_field(
    data.field[0].cutplane(1 * mm, x=0).scalar("Ey", "complex"), symmetrize=True
)
model.display.show()
