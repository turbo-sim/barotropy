import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

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
case_indices = [6]
case_table = case_table[case_table["index"].isin(case_indices)]

row = case_table.iloc[0]
efficiency_array = np.linspace(0.001, 1, 6)

# Define plot settings
save_figures = True
colors = mpl.cm.Blues(np.linspace(0.2, 0.8, len(efficiency_array)))
figs_axes = {}
plot_settings = [
    {
        "var_x": "pressure",
        "var_y": "density",
        # "y_lim": [1, 2000],
        # "y_scale": "log",
        "xlabel": "Pressure (bar)",
        "ylabel": "Density (kg/m$^3$)",
    },
    {
        "var_x": "pressure",
        "var_y": "temperature",
        "xlabel": "Pressure (bar)",
        "ylabel": "Temperature (K)",
    },
    {
        "var_x": "pressure",
        "var_y": "void_fraction",
        "y_lim": [0, 1],
        "xlabel": "Pressure (bar)",
        "ylabel": "Void Fraction (-)",
    },
    {
        "var_x": "pressure",
        "var_y": "speed_of_sound",
        "xlabel": "Pressure (bar)",
        "ylabel": "Speed of Sound (m/s)",
    },
]


# Loop over cases
for i, efficiency in enumerate(efficiency_array):

    # Process new row
    print(f"Processing case {i+1} of {len(efficiency_array)}")
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
        efficiency=efficiency,
        ODE_solver="LSODA",
        ODE_tolerance=1e-9,
        polynomial_degree=8,
        polynomial_format="horner",
        output_dir=dir_out,
    )

    # Evaluate barotropic model and export polynomial expressions
    model.solve()

    # Plot data for each variable
    for j, settings in enumerate(plot_settings):
        figs_axes = bpy.plot_xy_data(
            data=model.states,
            figs_axes=figs_axes,
            settings=settings,
            figsize=(6.0, 5.0),
            label=rf"$\eta={100*efficiency:0.0f}$ %",
            color=colors[i],
        )

# Adjust figures and create legends
for key, (fig, ax) in figs_axes.items():
    fig.tight_layout(pad=1)
    ax.legend(loc="best")

    if save_figures:
        file_path = os.path.join(DIR_OUT, f"sensitivity_efficiency_{key[0]}_{key[1]}")
        bpy.savefig_in_formats(fig, file_path)

# Show figure
plt.show()
