import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import barotropy as bpy

bpy.print_package_info()
bpy.set_plot_options()

# TODO: ideas to improve fitting. Doing fitting in log(density) scale
# TODO: double precission arithmetic in CFX
# TODO: custom fitting with non-negativity constraint using non-linear solver

from matplotlib import cm

# Create folder to save results
DIR_OUT = "./output"
# DIR_OUT = "./barotropic_model_20240805"
# DIR_OUT = "./barotropic_model_20240820"
if not os.path.isdir(DIR_OUT):
    os.makedirs(DIR_OUT)

# Load cases from Excel
case_table = pd.read_excel("./simulation_cases.xlsx")
case_indices = case_table["index"]
case_indices = [11, 12, 13]
# case_indices = [2]
case_table = case_table[case_table["index"].isin(case_indices)]

# Define plot settings
var_x = "s"
var_y = "T"
save_figures = True
colors = cm.magma(np.linspace(0.2, 0.8, max(len(case_indices) - 1 , 1)))  # Reds, Blues, magma

# Loop over cases
for i, (idx, row) in enumerate(case_table.iterrows()):

    # Process new row
    print(f"Processing case {i+1} of {len(case_table)}")
    print(row)
    print()
    color = "black" if i == 0 else colors[i - 1]
    print(row["tag"])
    dir_out = os.path.join(DIR_OUT, row["tag"])

    # Create barotropic model object
    fluid_name = row["fluid_name"]
    model = bpy.BarotropicModel(
        fluid_name=fluid_name,
        T_in=row["T0_in"],
        p_in=row["p0_in"],
        p_out=row["p_out"],
        efficiency=row["polytropic_efficiency"],
        calculation_type=row["calculation_type"],
        blending_onset=row["q_onset"],
        blending_width=row["q_transition"],
        # blending_width=0.05,
        HEOS_solver="hybr",
        ODE_solver="LSODA",
        ODE_tolerance=1e-9,
        polynomial_degree=8,
        polynomial_format="horner",
        output_dir=dir_out,
    )

    # Evaluate barotropic model and export polynomial expressions
    model.solve()
    model.fit_polynomials()
    model.export_expressions_fluent(output_dir=dir_out)
    model.export_expressions_cfx(output_dir=dir_out)
    for var in model.poly_fitter.variables:
        model.poly_fitter.plot_polynomial_and_error(var=var, savefig=True, showfig=True)

    # Create figure for plotting
    if i == 0:
        # Plot phase diagram
        fluid = bpy.Fluid(name=fluid_name, backend="HEOS", exceptions=True)
        fig, (ax_1, ax_2) = plt.subplots(
            1, 2, figsize=(12.0, 5.0), gridspec_kw={"wspace": 0.25}
        )
        fig.suptitle("Smooth barotropic model for expansion of nitrogen (LAES turbine case)", fontsize=14)
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
        color=color,
        linewidth=1.25,
        label=rf"$q_\text{{onset}}={row['q_onset']:0.2f}$, {row['calculation_type']}",
    )

    ax_2.set_yscale("log")
    ax_2.plot(
        model.states["p"],
        model.states["rho"],
        # linewidth=1.0,
        marker="o",
        markersize=3,
        color=color,
        label=rf"$q_\text{{onset}}={row['q_onset']:0.2f}$, {row['calculation_type']}",
    )

# Add legend
ax_1.legend(loc="lower right")
ax_2.legend(loc="lower right")
# bpy.scale_graphics_x(fig, +1e-3, mode="multiply")
# bpy.scale_graphics_y(fig, -273.15, mode="add")
# fig.tight_layout(pad=1)

# Save figures
if save_figures:
    file_path = os.path.join(DIR_OUT, f"barotropic_model_{var_x}_{var_y}_diagram")
    bpy.savefig_in_formats(fig, file_path)

# Show figure
plt.show()
