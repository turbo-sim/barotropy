
# %%

import os
import pandas as pd
import matplotlib.pyplot as plt

import barotropy as bpy


CASE_INDEX = 3

# Define plot settings
bpy.print_package_info()
bpy.set_plot_options()
SHOW_FIG = True
SAVE_FIG = True

# Create folder to save results
DIR_OUT = "./output"
save_figures = True
if not os.path.isdir(DIR_OUT):
    os.makedirs(DIR_OUT)

# Read the case summary file
EXCEL_DATA_SUMMARY = 'nakagawa_2009_cases.xlsx'
# EXCEL_DATA_PROFILES = 'pressure_profile_summary_converted.xlsx'
data_summary = pd.read_excel(EXCEL_DATA_SUMMARY)
# data_profiles = pd.read_excel(EXCEL_DATA_PROFILES)
case_data = data_summary[data_summary["Case"] == CASE_INDEX].squeeze()

# Create fluid object
fluid_name = data_summary['fluidname'].values[0]
fluid = bpy.Fluid(name=fluid_name, backend="HEOS")

# Define case parameters
dT_subcooling = 5

# total preessure in Pa
p_in = float(case_data["P_0_in"])

# Static pressure in Pa
p_out = case_data["P_out"]
# p_out = p_in/1.5

# Total temperature in K
T_in = case_data["T_0_in"] + 273.15

calculation_type = case_data["calculation_type"]
q_onset = case_data["q_onset"] # use the liquid line till quality of 5%
q_transition = case_data["q_transition"] # transient from metastable liquid to equilibrium zone
polytropic_efficiency = 1.00

print("Blending options:")
print(f'        calculation type: {calculation_type}')
print(f'        q_onset: {q_onset}')
print(f'        q_transition: {q_transition}')

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
    polynomial_degree=5,
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
    fluid=fluid, var_x="s", var_y="T", savefig=SHOW_FIG, showfig=SAVE_FIG
)
for var in model.poly_fitter.variables:
    model.poly_fitter.plot_polynomial_and_error(
        var=var, savefig=SHOW_FIG, showfig=SAVE_FIG
    )

# Show figures
plt.show()
