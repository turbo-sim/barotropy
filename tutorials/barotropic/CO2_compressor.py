import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# TODO, fix barotropy for the CO2 condensation cases. I am getting non-convergence
# TODO, fix rhe slope discontinuity, do the optimization of the polynomials myself


import barotropy as bpy

bpy.print_package_info()
bpy.set_plot_options()


# Create folder to save results
DIR_OUT = "./output_CO2"
save_figures = True
if not os.path.isdir(DIR_OUT):
    os.makedirs(DIR_OUT)

# Create fluid object
fluid_name = "CO2"
fluid = bpy.Fluid(name=fluid_name, backend="HEOS")

# Define case parameters
T_in = 32.47 + 273.15  # Degree Kelvin
p_in = 75.8*1e5        # Pa

state_in = fluid.get_state(bpy.PT_INPUTS, p_in, T_in)
s_in = state_in.s
state_in = fluid.get_state(bpy.PSmass_INPUTS, p_in*2, s_in)

p_in = state_in.p
T_in = state_in.T
p_out = fluid.critical_point.p/0.99


polytropic_efficiency = 1.00
calculation_type = "equilibrium"
q_onset = 0.99
q_transition = 0.02
# state_sat = fluid.get_state(bpy.PQ_INPUTS, p_in, 0)
# state_in = fluid.get_state(bpy.PT_INPUTS, state_sat.p, state_sat.T-dT_subcooling)


# Create barotropic model object
model = bpy.BarotropicModel(
    fluid_name=fluid_name,
    T_in=T_in,
    p_in=p_in,
    p_out=p_out,
    efficiency=polytropic_efficiency,
    calculation_type=calculation_type,
    blending_onset=q_onset,
    blending_width=q_transition,
    HEOS_solver="hybr",
    ODE_solver="LSODA",
    ODE_tolerance=1e-9,
    polynomial_degree=8,
    polynomial_format="horner",
    output_dir=DIR_OUT,
)

# Evaluate barotropic model and export polynomial expressions
model.solve()
model.fit_polynomials()
model.export_expressions_fluent(output_dir=DIR_OUT)
# model.export_expressions_cfx(output_dir=DIR_OUT)
# for var in model.poly_fitter.variables:
#     model.poly_fitter.plot_polynomial_and_error(var=var, savefig=True, showfig=True)


# Plot phase diagram
var_x = "s"
var_y = "T"
fig, (ax_1, ax_2) = plt.subplots(1, 2, figsize=(12.0, 5.0), gridspec_kw={"wspace": 0.25})
fig.suptitle(f"Barotropic model for expansion of {fluid_name}", fontsize=14)
ax_1.set_xlabel("Entropy (J/kg/K)\n")
ax_1.set_ylabel("Temperature (K)")
ax_2.set_xlabel("Pressure (Pa)\n")
ax_2.set_ylabel("Density (kg/m$^3$)")
# s_min, s_max = 1.3, 1.7
# T_min, T_max = 25, 35
# ax.set_xlim([s_min, s_max])
# ax.set_ylim([T_min, T_max])
ax_1 = fluid.plot_phase_diagram(
    var_x,
    var_y,
    axes=ax_1,
    plot_critical_point=True,
    plot_saturation_line=True,
    plot_spinodal_line=True,
    plot_quality_isolines=True,
    N=50,
)

ax_1.plot(
    model.states[var_x],
    model.states[var_y],
    # linewidth=1.0,
    # marker="o",
    # markersize=3,
    color=bpy.COLORS_MATLAB[0],
    linewidth=1.25,
)


ax_2.plot(
    model.states["p"],
    model.states["density"],
    linewidth=1.25,
    marker="none",
    markersize=3,
    linestyle="-",
    color=bpy.COLORS_MATLAB[0],
    label="Calculated data"
)
p = np.linspace(p_out, p_in, 1000)
# ax_2.plot(
#     p,
#     model.poly_fitter.evaluate_polynomial(p, "density"),
#     linewidth=1.25,
#     marker="none",
#     markersize=3,
#     linestyle="-.",
#     color=bpy.COLORS_MATLAB[1],
#     label="Polynomial fit"
# )

ax_2.legend(loc="lower right")


# Show figure
plt.show()
