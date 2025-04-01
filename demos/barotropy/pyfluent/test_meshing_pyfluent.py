import os
if not os.getenv('FLUENT_PROD_DIR'):
    import ansys.fluent.core as pyfluent
    flglobals = pyfluent.setup_for_fluent(product_version="24.2.0", mode="meshing", version="3d", precision="double", processor_count=10)
    globals().update(flglobals)

meshing.tui.file.import_.cad_geometry('yes', '"barana_geom.agdb"')
meshing.tui.file.import_.cad_geometry('yes', r'"C:\Users\alber\OneDrive - Danmarks Tekniske Universitet\Desktop\barana_geom.agdb"', 'mm', 'cad-faceting', 'no')
meshing.tui.file.import_.cad_geometry('yes', r'"C:\Users\alber\OneDrive - Danmarks Tekniske Universitet\Desktop\barana_geom.agdb"', 'mm', 'cad-faceting', 'no')
workflow.InitializeWorkflow(WorkflowType=r'2D Meshing')
meshing.tui.file.import_.cad('yes', r'"C:\Users\alber\OneDrive - Danmarks Tekniske Universitet\Desktop\barana_geom.agdb"', 'no', 'yes', '0', 'yes', 'mm', 'ok')
workflow.InitializeWorkflow(WorkflowType=r'2D Meshing')
meshing.tui.boundary.feature.list_edge_zones('()')
meshing.exit()
