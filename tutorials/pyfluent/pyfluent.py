import barotropy as bpy

import ansys.fluent.core as pyfluent
meshing = pyfluent.launch_fluent(mode=pyfluent.FluentMode.MESHING, product_version=pyfluent.FluentVersion.v242)
watertight = meshing.watertight()
watertight.import_geometry.file_name = pyfluent.examples.download_file("mixing_elbow.pmdb","pyfluent/mixing_elbow")
watertight.import_geometry()
watertight.create_volume_mesh()
meshing.switch_to_solver()



# setup, solution = solver.settings.setup, solver.settings.solution
# setup.boundary_conditions.set_zone_type(zone_list=["cold-inlet", "hot-inlet"], new_type="velocity-inlet")
# setup.boundary_conditions.set_zone_type(zone_list=["outlet"], new_type="pressure-outlet")
# setup.cell_zone_conditions.set_zone_type(zone_list="elbow-fluid", new_type="fluid")
# solution.initialization.hybrid_initialize()
# solution.run_calculation.iterate(iter_count=100)
# velocity_data = solver.fields.field_data.get_vector_field_data(field_name="velocity", surface_name="cold-inlet")