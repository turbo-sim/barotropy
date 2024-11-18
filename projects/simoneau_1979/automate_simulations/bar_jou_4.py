import os
if not os.getenv('FLUENT_PROD_DIR'):
    import ansys.fluent.core as pyfluent
    flglobals = pyfluent.setup_for_fluent(product_version="24.2.0", mode="solver", version="2d", precision="double", processor_count=10)
    globals().update(flglobals)

solver.tui.plot.plot('yes', '"test_plot_pressure.csv"', 'yes', 'no', 'no', 'pressure', 'no', 'no', 'x-coordinate', 'wall', '()')
