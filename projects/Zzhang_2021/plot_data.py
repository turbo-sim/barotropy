import os
import copy
from parse import parse
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import barotropy as bpy


# Define plot settings
save_figures = True
show_figures = True
case_name = "zhang_dykas_2021"
case_name_bis = "Zhang et al. (2021)"
bpy.set_plot_options(fontsize=14, grid=False)

# Check if the 'figures' directory exists. If not, create it.
figures_path = "output"
if not os.path.exists(figures_path):
    os.makedirs(figures_path)


# Read data from the excel file
nozzle_geom = pd.read_excel("validation_data.xlsx", sheet_name="nozzle_coordinates")
case_data = pd.read_excel("validation_data.xlsx", sheet_name="case_definition")
exp_data = pd.read_excel("validation_data.xlsx", sheet_name="experimental_data")

# Define the data of each case
fluid_name = "water"
fluid = bpy.Fluid(fluid_name)
T_in = case_data["T0_in[degC]"].to_numpy() + 273.15
p_in = case_data["p0_in[kPa]"].to_numpy() * 1e3
area_ratio = nozzle_geom["Y[m]"].to_numpy()[-1] / np.min(nozzle_geom["Y[m]"].to_numpy())

# Plot the nozzle coordinates
# plt.figure(figsize=(10, 6))
fig, ax = plt.subplots()
ax.set_xlabel("$x$ coordinates [mm]")
ax.set_ylabel("$y$ coordinates [mm]")
ax.set_title(f"{case_name_bis} nozzle coordinates")
ax.grid(visible=True)
ax.axis("equal")  # Ensure equal scaling on both axes
ax.plot(
    nozzle_geom["X[m]"] * 1e3,
    +nozzle_geom["Y[m]"] * 1e3,
    "-o",
    color="black",
    markersize=2.5,
    label="Nozzle Profile",
    linewidth=1.0,
    markeredgewidth=1.0
)
ax.plot(
    nozzle_geom["X[m]"] * 1e3,
    -nozzle_geom["Y[m]"] * 1e3,
    "-o",
    color="black",
    markersize=2.5,
    linewidth=1.0,
    markeredgewidth=1.0
)
ax.legend(loc="upper right")
fig.tight_layout(pad=1.00)

# Save plots
if save_figures:
    filepath = os.path.join(figures_path, f"{case_name}_nozzle_coordinates")
    bpy.savefig_in_formats(fig, filepath)

# Plot the nozzle expansions in a p-s diagram
fig, ax = plt.subplots()
ax.set_title(f"{case_name_bis} nozzle expansions")
ax.set_xlabel("$s$ - Entropy (kJ/kg)")
ax.set_ylabel("$p$ - Pressure (kPa)")
# ax.set_xlim([0.7, 1.1])
# ax.set_ylim([0, 2.2])
ax.set_yscale("log")
ax.grid(visible=False, which="both")

# Plot the phase diagram of the fluid
prop_x = "s"
prop_y = "p"
quality_isolines = np.linspace(0, 1, 21)
quality_delta = quality_isolines[1] - quality_isolines[0]
fluid.plot_phase_diagram(
    prop_x,
    prop_y,
    axes=ax,
    plot_critical_point=True,
    plot_saturation_line=True,
    plot_quality_isolines=True,
    quality_labels=False,
    quality_levels=quality_isolines
)

# Loop over all cases
for i, (T, p) in enumerate(zip(T_in, p_in), start=1):
    # Compute the exit pressure
    state_outlet, state_sonic, state_inlet = bpy.compute_supersonic_exit(
        T, p, area_ratio, fluid
    )

    # Plot the corresponding isentropic expansion
    x = [state_inlet[prop_x], state_outlet[prop_x]]
    y = [state_inlet[prop_y], state_outlet[prop_y]]
    ax.plot(
        x,
        y,
        "-o",
        markersize=3,
        linewidth=0.75,
        markeredgewidth=0.75,
        label=f"Case {i}",
    )

ax.plot([], [], color="black", linestyle=":", label=rf"$\Delta q={quality_delta*100}\%$")

# Re-scale the plot
ax.legend(loc="upper right", fontsize=10)
bpy.scale_graphics_x(fig, 1 / 1e3, mode="multiply")
bpy.scale_graphics_y(fig, 1 / 1e3, mode="multiply")
ax.relim()
fig.tight_layout(pad=1, w_pad=None, h_pad=None)

# Save plots
if save_figures:
    filepath = os.path.join(figures_path, f"{case_name}_ps_diagram")
    bpy.savefig_in_formats(fig, filepath)


# Plot the pressure distributions
fig = plt.figure(figsize=(6.4, 4.8))
ax = fig.gca()
ax.set_xlabel("Axial position (mm)")
ax.set_ylabel("Static pressure (kPa)")
ax.set_xscale("linear")
ax.set_yscale("linear")
for case_index in exp_data['Case'].unique():

    # Get the case specifications
    T0_in = case_data[case_data['Case'] == case_index]["T0_in[degC]"].iloc[0]
    p0_in = case_data[case_data['Case'] == case_index]["p0_in[kPa]"].iloc[0]

    # Plot experimental data
    _exp_data = exp_data[exp_data['Case'] == case_index]
    ax.plot(
        _exp_data["X[m]"]*1000,
        _exp_data["p[kPa]"],
        marker="o",
        linewidth=1.25,
        # color="red",
        linestyle="-",
        label=f"Case {case_index}",
    )

# Create legend
leg = ax.legend(loc="upper right")

# Adjust PAD
fig.tight_layout(pad=1, w_pad=None, h_pad=None)

# Save plots
if save_figures:
    filepath = os.path.join(figures_path, f"{case_name}_experimental_data")
    bpy.savefig_in_formats(fig, filepath)

# Display figures
if show_figures:
    plt.show()
