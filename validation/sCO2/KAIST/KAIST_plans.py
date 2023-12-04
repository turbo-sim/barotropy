import os
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd

import barotropy as bpy
import barotropy.fluid_properties as props
import barotropy.sCO2_utilities as sco2


# ---------------------------------------------------------------------------- #
# Load and postprocess data
# ---------------------------------------------------------------------------- #
OUTPUT_DIR = "results"
SAVE_FIGS = True
bpy.set_plot_options(minor_ticks=False)

# Create output directory if it does not exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Create fluid object
fluid = props.Fluid(name="CO2", backend="HEOS", exceptions=True)

# Planned experiments
p_in = np.asarray([75e5, 80e5])/ 1e5
d_in = fluid.critical_point.d*np.asarray([0.8, 0.9, 1.0, 1.1, 1.2])
combinations = list(itertools.product(p_in, d_in))
df = pd.DataFrame(combinations, columns=["pressure_in", "density_in"])
df = df.assign(**df.apply(lambda row: sco2.compute_inlet_state(row, fluid), axis=1))
df["compressor"] = df["pressure_in"].apply(lambda p: "Near saturation" if p == 75 else "High pressure")
df.to_excel("KAIST_plans.xlsx")
df["s_in"] *= 1e-3
df["p_in"] *= 1e-5
df["T_in"] -= 273.15

# Define plot settings
bpy.set_plot_options(minor_ticks=False, grid=False)
markers = ["P", "v", "s", "^", "o", "X", '<']
marker_sizes = {'o': 3.5, 'v': 4.5, 's': 3.5, '^': 4.5, 'P': 4.5, 'X': 4.5, '<': 3.5}
n = len(df["compressor"].unique())
colors = plt.get_cmap("magma")(np.linspace(0.80, 0.3, n))


# ---------------------------------------------------------------------------- #
# Plot inlet states on p-s diagram
# ---------------------------------------------------------------------------- #
# Create figure
fig, ax = plt.subplots(figsize=(6.8, 4.8))
ax.set_xlabel("Entropy [kJ/kg/K]")
ax.set_ylabel("Pressure [bar]")
range_x = np.linspace(1000, 2000, 50)
range_y = np.linspace(55, 105, 51)*1e5
range_z = np.asarray([30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45]) + 273.15
ax.set_xlim(np.asarray([range_x[0], range_x[-1]]) / 1e3)
ax.set_ylim(np.asarray([range_y[0], range_y[-1]]) / 1e5)

# Plot diagram
ax = sco2.plot_ps_diagram(fluid, ax, range_x, range_y, range_z)

# Loop over all cases
for j, compressor in enumerate(sorted(df["compressor"].unique())):
    subset = df[df["compressor"] == compressor]
    ax.plot(
        subset["s_in"] * 1e3,
        subset["p_in"] * 1e5,
        markerfacecolor="w",
        linewidth=0.75,
        linestyle="none",
        markeredgewidth=0.75,
        marker=markers[j % len(markers)],
        color=colors[j % len(colors)],
        markersize=marker_sizes[markers[j % len(markers)]],
        label=f"{compressor}",
    )

# Place the legend outside of the figure
ax.legend(loc="upper right", fontsize=10)
fig.tight_layout(pad=1)

# Scale to nice units
bpy.scale_graphics_x(fig, +1e-3, mode="multiply")
bpy.scale_graphics_y(fig, +1e-5, mode="multiply")

# Save figure
if SAVE_FIGS:
    base_filename = os.path.join(OUTPUT_DIR, f"KAIST_plans_data_ps_diagram")
    bpy.savefig_in_formats(fig, base_filename, formats=[".png", ".svg"])


# ---------------------------------------------------------------------------- #
# Plot inlet states on T-s diagram
# ---------------------------------------------------------------------------- #
# Create figure
fig, ax = plt.subplots(figsize=(6.8, 4.8))
ax.set_xlabel("Entropy [kJ/kg/K]")
ax.set_ylabel("Temperature [$^\circ$C]")
range_x = np.linspace(1000, 2000, 50)
range_y = np.linspace(5, 55, 50) + 273.15
range_z = np.asarray([50, 55, 60, 65, 70, 75, 80, 85, 90]) * 1e5
ax.set_xlim(np.asarray([range_x[0], range_x[-1]]) / 1e3)
ax.set_ylim(np.asarray([range_y[0], range_y[-1]]) - 273.15)

# Plot diagram
ax = sco2.plot_Ts_diagram(fluid, ax, range_x, range_y, range_z)

# Loop over all cases
for j, compressor in enumerate(sorted(df["compressor"].unique())):
    print(j, compressor)
    subset = df[df["compressor"] == compressor]
    ax.plot(
        subset["s_in"] * 1e3,
        subset["T_in"] + 273.15,
        markerfacecolor="w",
        linewidth=0.75,
        linestyle="none",
        markeredgewidth=0.75,
        marker=markers[j % len(markers)],
        color=colors[j % len(colors)],
        markersize=marker_sizes[markers[j % len(markers)]],
        label=f"{compressor}",
    )

# Place the legend outside of the figure
ax.legend(loc="lower right", fontsize=10)
fig.tight_layout(pad=1)

# Scale to nice units
bpy.scale_graphics_x(fig, +1e-3, mode="multiply")
bpy.scale_graphics_y(fig, -273.15, mode="add")

# Save figure
if SAVE_FIGS:
    base_filename = os.path.join(OUTPUT_DIR, f"KAIST_plans_Ts_diagram")
    bpy.savefig_in_formats(fig, base_filename, formats=[".png", ".svg"])


# Show figures
plt.show()