import os
if not os.getenv('FLUENT_PROD_DIR'):
    import ansys.fluent.core as pyfluent
    flglobals = pyfluent.setup_for_fluent(product_version="24.2.0", mode="solver", version="2d", precision="double", processor_count=20)
    globals().update(flglobals)

solver.solution.initialization.init_flow_statistics()
solver.solution.initialization.initialization_type = "standard"
solver.solution.initialization.standard_initialize()
solver.solution.initialization.initialization_type = "standard"
solver.solution.initialization.standard_initialize()
solver.solution.initialization.initialization_type = "hybrid"
solver.solution.initialization.hybrid_initialize()
solver.solution.run_calculation.iterate(iter_count = 50)
solver.tui.plot.plot('yes', '"test_plot.csv"', 'yes', 'no', 'no', 'pressure', 'no', 'no', 'x-coordinate', 'wall', '0', '()')
