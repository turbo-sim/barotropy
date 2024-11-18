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

# Barotropic figures settings
SHOW_FIG = True
SAVE_FIG = True

# Read Case Summary
CASES = ["design_inlet_condition"]
EXCEL_DATAFILE = "./combined_data.xlsx"  # Case summary file
data = pd.read_excel(EXCEL_DATAFILE)
case_data = data[data["tag"].isin(CASES)]
case_data = data.iloc[0, :]

# Import experimental data and convert to SI units
p_in = case_data["PT-109 (bar) mean"] * 1e5
p_out = case_data["PT-110 (bar) mean"] * 1e5
T_in = case_data["TE-109 (degC) mean"] + 273.15
T_out = case_data["TE-110 (degC) mean"] + 273.15
mass_flow_in = case_data["FT-101F (kg/min) mean"] / 60
mass_flow_out = case_data["FT-102F (kg/min) mean"] / 60
RPM = case_data["FU,1 (hz) mean"] * 60 

print(f"Inlet pressure: {p_in:0.3f} Pa")
print(f"Inlet temperature: {T_in:0.3f} K")
print(f"Inlet mass flow: {mass_flow_in:0.3f} kg/s")
print(f"Outlet pressure: {p_out:0.3f} Pa")
print(f"Angular speed: {RPM:0.3f} RPM")


# Create fluid object
fluid_name = "CO2"
fluid = bpy.Fluid(name=fluid_name, backend="HEOS")

# Compute properties as expansion from the outlet pressure and inlet entropy
state_in = fluid.get_state(bpy.PT_INPUTS, p_in, T_in)
state_out_s = fluid.get_state(bpy.PSmass_INPUTS, p_out, state_in.s)
p_in = state_out_s.p
T_in = state_out_s.T
p_out = 0.7*fluid.critical_point.p

# Additional parameters for the calculations
polytropic_efficiency = 1.00
calculation_type = "equilibrium"
q_onset = 1.00
q_transition = 0.05

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
    polynomial_degree=3,
    polynomial_format="horner",
    output_dir=DIR_OUT,
    polynomial_variables=["density", "viscosity", "speed_sound"],
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
