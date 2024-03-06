import os
import copy
from parse import parse
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import barotropy as brtp


# Define plot settings
save_figures = True
show_figures = True
case_name = "zhang_dykas_2021"
case_name_bis = "Zhang et al. (2021)"
brtp.set_plot_options(fontsize=14, grid=False)

# Check if the 'figures' directory exists. If not, create it.
figures_path = "figures"
if not os.path.exists(figures_path):
    os.makedirs(figures_path)


# Read data from the excel file
nozzle_geom = pd.read_excel("validation_data.xlsx", sheet_name="nozzle_coordinates")
case_data = pd.read_excel("validation_data.xlsx", sheet_name="case_definition")
exp_data = pd.read_excel("validation_data.xlsx", sheet_name="experimental_data")

# Define the data of each case
fluid_name = "water"
fluid = brtp.Fluid(fluid_name)
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

if save_figures:
    plt.savefig(os.path.join(figures_path, f"{case_name}_nozzle_coordinates.png"), dpi=500)
    plt.savefig(os.path.join(figures_path, f"{case_name}_nozzle_coordinates.svg"))


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
brtp.plot_phase_diagram(
    prop_x,
    prop_y,
    fluid,
    plot_critical_point=True,
    plot_saturation_line=True,
    plot_triple_point=False,
    plot_quality_isolines=True,
    quality_labels=False,
    quality_levels=quality_isolines
)

# Loop over all cases
for i, (T, p) in enumerate(zip(T_in, p_in), start=1):
    # Compute the exit pressure
    state_outlet, state_sonic, state_inlet = brtp.compute_supersonic_exit(
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

ax.plot([], [], color="black", linestyle=":", label=f"$\Delta q={quality_delta*100}\%$")

# Re-scale the plot
ax.legend(loc="upper right", fontsize=10)
brtp.scale_graphics_x(fig, 1 / 1e3, mode="multiply")
brtp.scale_graphics_y(fig, 1 / 1e3, mode="multiply")
ax.relim()
fig.tight_layout(pad=1, w_pad=None, h_pad=None)

if save_figures:
    plt.savefig(os.path.join(figures_path, f"{case_name}_ps_diagram.png"), dpi=500)
    plt.savefig(os.path.join(figures_path, f"{case_name}_ps_diagram.svg"))


# Plot the pressure distributions
for case_index in exp_data['Case'].unique():

    # Get the case specifications
    T0_in = case_data[case_data['Case'] == case_index]["T0_in[degC]"].iloc[0]
    p0_in = case_data[case_data['Case'] == case_index]["p0_in[kPa]"].iloc[0]

    # Create figure
    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.gca()
    ax.set_title(f"{case_name_bis} - Case {case_index}")
    ax.set_xlabel("Axial position (mm)")
    ax.set_ylabel("Static pressure (kPa)")
    ax.set_xscale("linear")
    ax.set_yscale("linear")
    # ax.set_ylim([0, np.ceil(pressure / 10) * 10])
    # ax.set_ylim([])
    # ax.set_xticks([])
    # ax.set_yticks([])

# # Plot simulation data
#     # Load Fluent XY datafile
# file_sim = os.path.join(results_dir, f"case_{case_index}_pressure.xy")
# df_sim = parse_fluent_xy(file_sim)
# order = df_sim["wall"]["Position"].argsort()
# x = df_sim["wall"]["Position"][order] * 1e3
# y = df_sim["wall"]["Static Pressure"][order] / 1e5
# ax.plot(
#     x,
#     y,
#     # color="black",
#     linewidth=1.5,
#     linestyle="-",
#     marker="none",
#     label="Simulation results",
# )

    # Plot experimental data
    _exp_data = exp_data[exp_data['Case'] == case_index]
    ax.plot(
        _exp_data["X[m]"]*1000,
        _exp_data["p[kPa]"],
        marker="o",
        linewidth=1.25,
        # color="red",
        linestyle="None",
        label="Experimental data",
    )

    # Create legend
    leg = ax.legend(loc="upper right")

    # Adjust PAD
    fig.tight_layout(pad=1, w_pad=None, h_pad=None)

    # Save plots
    if save_figures:
        plt.savefig(os.path.join(figures_path, f"{case_name}_validation_case_{case_index}.png"), dpi=500)
        plt.savefig(os.path.join(figures_path, f"{case_name}_validation_case_{case_index}.svg"))



# Display figures
if show_figures:
    plt.show()
