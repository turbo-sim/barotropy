import os
import pandas as pd
import matplotlib.pyplot as plt

import barotropy as bpy
bpy.print_package_info()
bpy.set_plot_options(fontsize=14, grid=False)

# Create folder to save results
DIR_OUT = "./unsteady_paper"
if not os.path.isdir(DIR_OUT):
    os.makedirs(DIR_OUT)

# Load cases from Excel
case_table = pd.read_excel("./simulation_cases.xlsx")
case_indices = case_table["index"]
case_indices = [6]
case_table = case_table[case_table["index"].isin(case_indices)]
for i, (idx, row) in enumerate(case_table.iterrows()):

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
    )

# Evaluate barotropic model and export polynomial expressions
model.solve()

# Create a figure and axis
fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(10, 5))
pad_x = 5
pad_y = 8

# Density
ax = axes[0][0]
ax.set_xlabel(r"Pressure (kPa)", labelpad=pad_x)
ax.set_ylabel(r"Density (kg/m$^3$)", labelpad=pad_y)
ax.set_yscale("log")
ax.plot(model.states["p"]/1e3, model.states["rho"], color="black")
ax.plot(model.states["p"]/1e3, model.states["rho_1"], color=bpy.COLORS_MATLAB[0])
ax.plot(model.states["p"]/1e3, model.states["rho_2"], color=bpy.COLORS_MATLAB[1])

# Void fraction
ax = axes[0][1]
ax.set_xlabel(r"Pressure (kPa)", labelpad=pad_x)
ax.set_ylabel(r"Void fraction (%)", labelpad=pad_y)
ax.plot(model.states["p"]/1e3, 100*model.states["void_fraction"], color="black")
ax.plot([], [], color="black", label="Mixture value")
ax.plot([], [], color=bpy.COLORS_MATLAB[0], label="Water value")
ax.plot([], [], color=bpy.COLORS_MATLAB[1], label="Nitrogen value")
ax.legend(loc="upper right", fontsize=11)

# Temperature
ax = axes[1][0]
ax.set_xlabel(r"Pressure (kPa)", labelpad=pad_x)
ax.set_ylabel(r"Temperature ($^\circ$C)", labelpad=pad_y)
ax.plot(model.states["p"]/1e3, model.states["T"] - 273.15, color="black")

# Speed of sound
ax = axes[1][1]
ax.set_xlabel(r"Pressure (kPa)", labelpad=pad_x)
ax.set_ylabel(r"Speed of sound (m/s)", labelpad=pad_y)
ax.set_yscale("log")
ax.plot(model.states["p"]/1e3, model.states["speed_sound"], color="black")
ax.plot(model.states["p"]/1e3, model.states["speed_sound_1"], color=bpy.COLORS_MATLAB[0])
ax.plot(model.states["p"]/1e3, model.states["speed_sound_2"], color=bpy.COLORS_MATLAB[1])

# Save figure
fig.tight_layout(pad=1.0)
filename = os.path.join(DIR_OUT, "barotropic_model")
bpy.savefig_in_formats(fig, filename)

# Show figure
plt.show()

