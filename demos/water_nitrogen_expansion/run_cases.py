import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import barotropy as bpy
bpy.print_package_info()
bpy.set_plot_options()

# Create folder to save results
DIR_OUT = "./output"
if not os.path.isdir(DIR_OUT):
    os.makedirs(DIR_OUT)

# Load cases from Excel
case_table = pd.read_excel("./simulation_cases.xlsx")
case_indices = case_table["index"]
# case_indices = [1]
case_table = case_table[case_table["index"].isin(case_indices)]

# Define plot settings
save_figures = True
colors = mpl.cm.Blues(np.linspace(0.2, 0.8, len(case_indices)))  # Reds, Blues, magma
figs_axes = {}
plot_settings = [
    {
        "var_x": "p",
        "var_y": "density",
        # "x_lim": [0, 20e5],
        "y_lim": [1, 2000],
        "y_scale": "log",
        "xlabel": "Pressure (bar)",
        "ylabel": "Density (kg/m$^3$)"
    },
    {
        "var_x": "p",
        "var_y": "T",
        "x_lim": [0, 20e5],
        "xlabel": "Pressure (bar)",
        "ylabel": "Temperature (K)"
    },
    {
        "var_x": "p",
        "var_y": "void_fraction",
        "x_lim": [0, 20e5],
        "y_lim": [0, 1],
        "xlabel": "Pressure (bar)",
        "ylabel": "Void Fraction (-)"
    },
    {
        "var_x": "p",
        "var_y": "speed_sound",
        "x_lim": [0, 20e5],
        "xlabel": "Pressure (bar)",
        "ylabel": "Speed of Sound (m/s)"
    }
]



# Loop over cases
for i, (idx, row) in enumerate(case_table.iterrows()):

    # Process new row
    print(f"Processing case {i+1} of {len(case_table)}")
    print(row)
    print()
    dir_out = os.path.join(DIR_OUT, f"case_{row['index']}")

    # Create barotropic model object
    model = bpy.BarotropicModel(
        fluid_name=[row["fluid_1"], row["fluid_2"]],
        mixture_ratio=row["mixture_ratio"],
        T_in=row["T0_in"],
        p_in=row["p0_in"],
        p_out=row["p_out"],
        efficiency=row["efficiency"],
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
        model.poly_fitter.plot_polynomial_and_error(var=var, savefig=True, showfig=False)

    # Plot data for each variable
    for j, settings in enumerate(plot_settings):
        figs_axes = bpy.plot_xy_data(
            data=model.states,
            figs_axes=figs_axes,
            settings=settings,
            figsize=(6.0, 5.0),
            label=rf"$R={row['mixture_ratio']:0.2f}$",
            color=colors[i],
        )
    
    # # Plot density for individual phases
    # if i == len(case_indices) - 1:
    #     fig, ax = figs_axes[('p', 'density')]
    #     ax.plot(model.states["p"], model.states["rho_1"], color="black", label="Water")
    #     ax.plot(model.states["p"], model.states["rho_2"], color=bpy.COLORS_MATLAB[1], label="Nitrogen")


# Adjust figures and create legends
for key, (fig, ax) in figs_axes.items():
    fig.tight_layout(pad=1)
    ax.legend(loc="best")

    if save_figures:
        file_path = os.path.join(DIR_OUT, f"barotropic_model_{key[0]}_{key[1]}")
        bpy.savefig_in_formats(fig, file_path)

# Show figures
if not os.environ.get("DISABLE_PLOTS"):
    plt.show()

