import os
import gmsh
import yaml
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import barotropy as bpy
import ansys.fluent.core as pyfluent

# Work done and to do:
#  - I created a class to generate nozzle meshes with Gmsh
#  - Functionality
#    - Create meshes with symmetry plane, mirrored nozzle or asymetric nozzle
#    - Create nozzle with and without inflation layers
#    - Complete control over the first cell element size and total number of cells (and approximate control over inflation ratio)
#    - Custom function to calculate and optimize the 1D inflation layers
#    - The class cna visualize and export the meshes
# - The .msh format from Gmsh is different and has nothing to do with the Fluent .msh format
# - The meshes exported in NASTRAN (.bdf) and ABAQUS (.inf) formats can be loaded in Fluent, but they lose the information about physical groups
# - It might be needed to postprocess the mesh in Fluent to rename the physical groups to assign boundary conditions
# - This is tricky because Fluent automatically names boundaries as "wall-1", "wall-2" etc upon loading and it is not trivial to know which edge is which
# - Using meshio to conver between mesh formats does not seem to help. Perhaps exporting to .cgns format could do the trick, but I have not managed to import .cgns mesh in fluent properly
# - One potential solution could be to import a NASTRAN or ABAQUS mesh in Fluent and then identify the coordinates of each element (e.g., its centroid) to be able to map them to the centroids of each element calculated in Gmsh
# - In this way we could change the name of the Fluent Zones and assign the proper boundary conditions automatically using PyFluent
# - This would provide a streamlined approach to automatically define geometry, boundary conditions and solver settings and have a simulation ready to be run in Ansys
# - The approach could easily be extended to blades or 3D geometries
# - I am not sure if the effort is worth it or if we should try to parameterize another meshing tool (e.g., Ansys mesher) to streamline the interface to Fluent/CFX


# Load the configuration file
with open("mesh_config.yaml", "r") as f:
    config = yaml.safe_load(f)

# Generate the mesh
out_filename = "nozzle.bdf"
mesher = bpy.NozzleMesher(config)
mesher.create_mesh()
mesher.show_mesh()
mesher.write_mesh(filename=out_filename)

# Launch Fluent
print("Launching Fluent...")
solver = pyfluent.launch_fluent(
    product_version=pyfluent.FluentVersion.v242,
    mode=pyfluent.FluentMode.SOLVER,
    ui_mode=pyfluent.UIMode.NO_GUI,
    dimension=pyfluent.Dimension.TWO,
    precision=pyfluent.Precision.DOUBLE,
    processor_count=1,
)
print("Fluent launched successfully.")

# Import mesh file
solver.file.import_.read(file_type="nastran-bulkdata", file_name=out_filename)

# solver.mesh.mesh_info(print_level=0)
# solver.mesh.modify_zones.list_zones()


df_zones = bpy.get_fluent_zone_table(solver)

print("Dataframe with zones info:")
print(df_zones)


