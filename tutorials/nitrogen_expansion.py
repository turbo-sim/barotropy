import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import barotropy as bpy

bpy.print_package_info()
bpy.set_plot_options()


# Create folder to save results
DIR_OUT = "./output"
save_figures = True
if not os.path.isdir(DIR_OUT):
    os.makedirs(DIR_OUT)


###################################################################################################################
# Select case number
cas = 69
###################################################################################################################


# Read Case Summary
case_summ = "../projects/Simoneau_Hendricks_1979/cases_summary.xlsx"
data = pd.read_excel(case_summ)
data = pd.DataFrame(data)
case_data = data[data.iloc[:, 0] == cas]

# Create fluid object
fluid_name = case_data["fluidname"].iloc[0] 
fluid = bpy.Fluid(name=fluid_name, backend="HEOS")


# Define case parameters
dT_subcooling = 10

# total preessure in Pa
p_in = float(case_data["P_0_in"])
# Static pressure in Pa
p_out = float(case_data["P_out"])
# Total temperature in K
T_in = float(case_data["T_0_in"])

polytropic_efficiency = 1.00
calculation_type = "blending"
q_onset = 0.10 # use the liquid line till quality of 5%
q_transition = 0.05 # transient from metastable liquid to equilibrium zone


# Create barotropic model object
# state_sat = fluid.get_state(bpy.PQ_INPUTS, p_in, 0)
# state_in = fluid.get_state(bpy.PT_INPUTS, state_sat.p, state_sat.T-dT_subcooling)
model = bpy.BarotropicModel(
    fluid_name=fluid_name,
    T_in=T_in,
    p_in=p_in,
    p_out=p_out/2,
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
    polynomial_variables=['density', 'viscosity', 'speed_sound']
)

# Evaluate barotropic model and export polynomial expressions
model.solve()
model.fit_polynomials()
model.export_expressions_fluent(output_dir=DIR_OUT)
model.export_expressions_cfx(output_dir=DIR_OUT)
for var in model.poly_fitter.variables:
    model.poly_fitter.plot_polynomial_and_error(var=var, savefig=True, showfig=True)



# print("hello world")

# expression = bpy.read_expression_file(filename="output/fluent_expressions.txt", vars=["density", "viscosity"])

# bpy.print_dict(expression)

# bpy.print_dict(expression)

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
ax_2.plot(
    p,
    model.poly_fitter.evaluate_polynomial(p, "density"),
    linewidth=1.25,
    marker="none",
    markersize=3,
    linestyle="-.",
    color=bpy.COLORS_MATLAB[1],
    label="Polynomial fit"
)

ax_2.legend(loc="lower right")


# Show figure
plt.show()
