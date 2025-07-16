import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import barotropy as bpy

bpy.print_package_info()
bpy.set_plot_options(grid=False)


# Create folder to save results
DIR_OUT = "./output"
save_figures = True
if not os.path.isdir(DIR_OUT):
    os.makedirs(DIR_OUT)

# Barotropic figures settings
SHOW_FIG = True
SAVE_FIG = False

# Create fluid object
fluid_name = "CO2"
fluid = bpy.Fluid(name=fluid_name, backend="HEOS")

# Define inlet state and outlet pressure
p_in = 2.0*fluid.critical_point.p
# s_in = 0.9*fluid.critical_point.s
s_in = 1.1*fluid.critical_point.s
state_in = fluid.get_state(bpy.PSmass_INPUTS, p_in, s_in)
p_out = 1.01*fluid.triple_point_liquid.p

# Additional parameters for the calculations
polytropic_efficiency = 1.0
q_onset = 0.95
q_transition = 0.05

# Create barotropic model object
model = bpy.BarotropicModel(
    fluid_name=fluid_name,
    T_in=state_in.T,
    p_in=state_in.p,
    p_min=p_out,
    p_max=p_in,
    efficiency=polytropic_efficiency,
    process_type="expansion",
    calculation_type="blending",
    blending_onset=q_onset,
    blending_width=q_transition,
    HEOS_solver="hybr",
    ODE_solver="LSODA",
    ODE_tolerance=1e-9,
    polynomial_degree=6,
    polynomial_format="horner",
    output_dir=DIR_OUT,
)

# Evaluate barotropic model and export polynomial expressions
model.solve()
model.fit_polynomials()
model.export_expressions_fluent(output_dir=DIR_OUT)
model.export_expressions_cfx(output_dir=DIR_OUT)
model.poly_fitter.plot_phase_diagram(
    fluid=fluid,
    var_x="s",
    var_y="T",
    savefig=SAVE_FIG,
    showfig=SHOW_FIG,
)
for var in model.poly_fitter.variables:
    model.poly_fitter.plot_polynomial_and_error(
        var=var,
        savefig=SAVE_FIG,
        showfig=SHOW_FIG,
    )


# Show figure
plt.show()
