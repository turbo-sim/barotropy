import os
import pandas as pd
import matplotlib.pyplot as plt

import barotropy as bpy

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

# Create fluid object
fluid_name = "nitrogen"
fluid = bpy.Fluid(name=fluid_name, backend="HEOS")

# # Read data from table
# case_index = 69
# case_summ = "cases_summary.xlsx"
# case_table = pd.read_excel(case_summ)
# case_data = case_table[case_table["Case"] == case_index].squeeze()
# T_in = case_data["T_0_in"]
# p_in = case_data["P_0_in"]
# p_out = case_data["P_out"]

# Define thermodynamic data manuallt
dT_subcooling = 4
p_in = 0.4 * fluid.critical_point.p
state_sat = fluid.get_state(bpy.PQ_INPUTS, p_in, 0)
state_in = fluid.get_state(bpy.PT_INPUTS, state_sat.p, state_sat.T - dT_subcooling)
T_in = state_in.T
p_out = fluid.triple_point_liquid.p

# Define barotropic model parameters
polytropic_efficiency = 1.00
calculation_type = "equilibrium"
q_onset = 0.05  # use the liquid line till quality of 5%
q_transition = 0.05  # transient from metastable liquid to equilibrium zone

# Create barotropic model object
model = bpy.BarotropicModel(
    fluid_name=fluid_name,
    T_in=T_in,
    p_in=p_in,
    p_out=p_out,
    efficiency=polytropic_efficiency,
    calculation_type=calculation_type,
    # blending_onset=q_onset,
    # blending_width=q_transition,
    HEOS_solver="hybr",
    ODE_solver="LSODA",
    ODE_tolerance=1e-9,
    polynomial_degree=8,
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
if not os.environ.get("DISABLE_PLOTS"):
    plt.show()
