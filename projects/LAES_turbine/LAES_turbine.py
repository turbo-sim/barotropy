import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import barotropy as bpy

bpy.print_package_info()
bpy.set_plot_options()


from matplotlib import cm

# Create folder to save results
DIR_OUT = "./output"
if not os.path.isdir(DIR_OUT):
    os.makedirs(DIR_OUT)

# Load cases from Excel
case_table = pd.read_excel("./simulation_cases.xlsx")
case_indices = case_table["index"]
# case_indices = [1, 2, 3]
case_indices = [2]
case_table = case_table[case_table["index"].isin(case_indices)]

# Define plot settings
var_x = "s"
var_y = "T"
n_points = None
# colors = bpy.COLORS_MATLAB
# colors = cm.Reds(np.linspace(0.5, 1.0, len(case_indices)))
colors = cm.magma(np.linspace(0.2, 0.8, len(case_indices)))

save_figures = True

# Loop over cases
for i, (idx, row) in enumerate(case_table.iterrows()):

    # Create fluid object
    fluid_name = row["fluid_name"]
    print(f"Processing case {i+1} of {len(case_table)}")
    print(row)
    print()
    fluid = bpy.Fluid(name=fluid_name, backend="HEOS", exceptions=True)

    p_in = row["p0_in"]
    p_out = row["p_out"]

    # Evaluate barotropic model
    polyfit = bpy.barotropic_model_one_fluid(
        fluid_name=fluid_name,
        T_in=row["T0_in"],
        p_in=p_in,
        p_out=p_out,
        efficiency=row["polytropic_efficiency"],
        calculation_type="blending",
        blending_onset=row["q_onset"],
        blending_width=row["q_transition"],
        n_eval=n_points,
        HEOS_solver="hybr",
        ODE_solver="BDF",
        ODE_tolerance=1e-9,
        polynomial_order=8,
        polynomial_vars=["density", "viscosity", "speed_sound"],
    )

    # # Create figure for plotting
    # if i == 0:
    #     # Plot phase diagram
    #     fig_1, ax_1 = plt.subplots(figsize=(6.0, 5.0))
    #     ax_1.set_xlabel("Entropy (J/kg/K)")
    #     ax_1.set_ylabel("Temperature (K)")
    #     # s_min, s_max = 1.3, 1.7
    #     # T_min, T_max = 25, 35
    #     # ax.set_xlim([s_min, s_max])
    #     # ax.set_ylim([T_min, T_max])
    #     ax_1 = fluid.plot_phase_diagram(
    #         var_x,
    #         var_y,
    #         axes=ax_1,
    #         plot_critical_point=True,
    #         plot_saturation_line=True,
    #         plot_spinodal_line=False,
    #         plot_quality_isolines=True,
    #     )

    # if i == 0:
    #     color = "black"
    # else:
    #     color = colors[i]
        
    # ax_1.plot(
    #     states[var_x],
    #     states[var_y],
    #     linewidth=1.0,
    #     # marker="o",
    #     # markersize=3,
    #     color=color,
    #     label=rf"$q_\text{{onset}}={row['q_onset']:0.2f}$, {row['multiphase_model']}",
    # )


    fig_2 = polyfit.plot_polynomial_and_error(var="density")
    # fig_3 = poly.plot_polynomial_and_error(var="viscosity")
    # fig_4 = poly.plot_polynomial_and_error(var="speed_sound")


    polyfit.export_fluent_expressions()



# # Adjust plot
# ax_1.legend(loc="lower right")
# ax_2.legend(loc="lower right")
# # bpy.scale_graphics_x(fig, +1e-3, mode="multiply")
# # bpy.scale_graphics_y(fig, -273.15, mode="add")
# fig_1.tight_layout(pad=1)
# fig_2.tight_layout(pad=1)


# # export_fluent_expressions(barotropic_model, fluid_name=fluid_nane, case_name=case_name, output_dir=output_folder)
# # export_cfx_expressions(barotropic_model, fluid_name=fluid_nane, case_name=case_name, output_dir=output_folder)


# # Save figures
# if save_figures:
#     bpy.savefig_in_formats(
#         fig_1, os.path.join(DIR_OUT, f"barotropic_model_{var_x}_{var_y}_diagram")
#     )
#     bpy.savefig_in_formats(fig_2, os.path.join(DIR_OUT, "barotropic_model_rho_p"))

# Show figure
plt.show()
