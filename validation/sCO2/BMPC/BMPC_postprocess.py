import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd

import barotropy as bpy
import barotropy.fluids.core as props
import barotropy.sCO2_utilities as sco2


# ---------------------------------------------------------------------------- #
# Load and postprocess data
# ---------------------------------------------------------------------------- #

# Postprocess raw data to:
#  1. Compute inlet state
#  2. Compute outlet isentropic state
#  3. Compute specific work
#  4. Compute actual mass flow rate from equivalent mass flow (when relevant)

# Define case parameters
CASE = 'BMPC'
SAVE_FIGS = True
OUTPUT_DIR = "results"
bpy.set_plot_options(minor_ticks=False)

# Create output directory if it does not exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Import data about simulation cases
df = pd.read_excel(
    f"{CASE}_data_digitized.xlsx",
    sheet_name="data",
    skiprows=lambda x: x in [1],  # Skip unit row
)

# Extract the desired dataset
df = df[(df["compressor"] == CASE) & (df["data_source"] == "EXP")]

# Create fluid
fluid = props.Fluid(name="CO2", backend="HEOS", exceptions=True)

# Compute reference state
p_ref = 91.7e5
T_ref = 35.556 + 273.15
reference_state = fluid.set_state(props.PT_INPUTS, p_ref, T_ref)

# Postprocess data
df = df.assign(**df.apply(lambda row: sco2.compute_inlet_state(row, fluid), axis=1))
df = df.assign(**df.apply(lambda row: sco2.compute_outlet_state_isentropic(row, fluid), axis=1))

# Compute 'mass_flow' for the datasets where 'mass_flow_eq' is provided
mask_actual = df['mass_flow'].isna() & df['mass_flow_eq'].notna()
df_actual = df[mask_actual].apply(lambda row: sco2.compute_actual_from_equivalent(row, reference_state, method="perfect_gas"), axis=1)
df.loc[mask_actual, df_actual.columns] = df_actual

# Compute 'mass_flow_eq' for the datasets where 'mass_flow' is provided
mask_eq = df['mass_flow_eq'].isna() & df['mass_flow'].notna()
df_eq = df[mask_eq].apply(lambda row: sco2.compute_equivalent_from_actual(row, reference_state, method="perfect_gas"), axis=1)
df.loc[mask_eq, df_eq.columns] = df_eq

# # Compute all the equivalent variables
# df = df.assign(**df.apply(lambda row: sco2.compute_equivalent_from_actual(row, reference_state, method="perfect_gas"), axis=1))

# Save postprocessed data
df.to_excel(f"{CASE}_data_postprocessed.xlsx")

# Compute equivalent variables
df_eq = df.assign(**df.apply(lambda row: sco2.compute_equivalent_from_actual(row, reference_state, method="real_gas"), axis=1))

# Define colors and markers for plots
markers = ["o", "^", "s", "v", "p", "+", "x"]
n = len(df["tag"].unique())
colors = plt.get_cmap("magma")(np.linspace(0.1, 0.8, n + 1))

# ---------------------------------------------------------------------------- #
# Plot raw performance data
# ---------------------------------------------------------------------------- #
# Create a figure and two subplots sharing the same legend
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 4.8), sharey=True)
fig.suptitle("sCO$_2$-HeRo compressor performance data", fontsize=14)
ax1.set_xlabel("Mass flow rate (kg/s)")
ax1.set_ylabel("Pressure ratio")
ax2.set_xlabel("Angular Speed (RPM)")
ax2.set_ylabel("Pressure ratio")
# ax1.set_xlim([0.0, 1.2])
# ax2.set_xlim([5000, 65000])
# ax1.set_ylim([0.95, 1.55])

# Loop over the unique values of 'tag'.
for i, tag in enumerate(sorted(df["tag"].unique())):
    # Extract the subset of data for each tag.
    subset = df[df["tag"] == tag]

    # Define the label
    label = f"Dataset {tag}"
    
    # Plot the subset on the first subplot
    ax1.plot(
        subset["mass_flow"],
        subset["pressure_ratio"],
        marker=markers[i % len(markers)],
        markersize=4.5,
        markerfacecolor="w",
        linewidth=0.75,
        linestyle="none",
        markeredgewidth=1.00,
        color=colors[i % len(colors)],
        label=label,
    )
    
    # Plot the subset on the second subplot
    ax2.plot(
        subset["angular_speed"],
        subset["pressure_ratio"],
        marker=markers[i % len(markers)],
        markersize=4.5,
        markerfacecolor="w",
        linewidth=0.75,
        linestyle="none",
        markeredgewidth=1.00,
        color=colors[i % len(colors)],
    )

# Adjust the layout and add legend
plt.subplots_adjust(right=0.75)
leg = fig.legend(loc='center left', bbox_to_anchor=(0.775, 0.5))
plt.subplots_adjust(bottom=0.2)

# Save figure
if SAVE_FIGS:
    base_filename = os.path.join(OUTPUT_DIR, f"{CASE}_performance_data")
    bpy.savefig_in_formats(fig, base_filename, formats=[".png", ".svg"])


# ---------------------------------------------------------------------------- #
# Plot equivalent performance data
# ---------------------------------------------------------------------------- #
# Create a figure and two subplots sharing the same legend
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 4.8), sharey=True)
fig.suptitle("sCO$_2$-HeRo compressor equivalent performance data", fontsize=14)
ax1.set_xlabel("Equivalent mass flow rate (kg/s)")
ax1.set_ylabel("Equivalent isentropic work (kJ/kg)")
ax2.set_xlabel("Equivalent angular Speed (RPM)")
ax2.set_ylabel("Equivalent isentropic work (kJ/kg)")
# ax1.set_xlim([0.0, 1.2])
# ax2.set_xlim([5000, 65000])

# Loop over the unique values of 'tag'.
for i, tag in enumerate(sorted(df_eq["tag"].unique())):
    # Extract the subset of data for each tag.
    subset = df_eq[df_eq["tag"] == tag]

    # Define the label
    label = f"Dataset {tag}"
    
    # Plot the subset on the first subplot
    ax1.plot(
        subset["mass_flow_eq"],
        subset["isentropic_work_eq"]/1e3,
        marker=markers[i % len(markers)],
        markersize=4.5,
        markerfacecolor="w",
        linewidth=0.75,
        linestyle="none",
        markeredgewidth=1.00,
        color=colors[i % len(colors)],
        label=label,
    )
    
    # Plot the subset on the second subplot
    ax2.plot(
        subset["angular_speed_eq"],
        subset["isentropic_work_eq"]/1e3,
        marker=markers[i % len(markers)],
        markersize=4.5,
        markerfacecolor="w",
        linewidth=0.75,
        linestyle="none",
        markeredgewidth=1.00,
        color=colors[i % len(colors)],
    )

# Adjust the layout and add legend
plt.subplots_adjust(right=0.75)
leg = fig.legend(loc='center left', bbox_to_anchor=(0.775, 0.5))
plt.subplots_adjust(bottom=0.2)

# Save figure
if SAVE_FIGS:
    base_filename = os.path.join(OUTPUT_DIR, f"{CASE}_performance_data_equivalent")
    bpy.savefig_in_formats(fig, base_filename, formats=[".png", ".svg"])
    df_eq.to_excel(f"{base_filename}.xlsx")



# ---------------------------------------------------------------------------- #
# Plot equivalent performance data
# ---------------------------------------------------------------------------- #
# Create a 3D plot
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111, projection='3d')
fig.suptitle("sCO$_2$-HeRo compressor equivalent performance data", fontsize=14)
ax.set_xlabel("Equivalent mass flow rate (kg/s)")
ax.set_zlabel("Equivalent isentropic work (kJ/kg)")
ax.set_ylabel("Equivalent angular speed (RPM)")

# Loop over the unique values of 'tag'.
for i, tag in enumerate(sorted(df_eq["tag"].unique())):
    # Extract the subset of data for each tag.
    subset = df_eq[df_eq["tag"] == tag]

    # Scatter plot for each subset
    ax.scatter(
        subset["mass_flow_eq"],
        subset["angular_speed_eq"],
        subset["isentropic_work_eq"]/1e3,
        marker=markers[i % len(markers)],
        edgecolors=colors[i % len(colors)],
        label=f"Dataset {tag}",
        facecolors="white",
        alpha=1.00,
    )

# Add a legend outside of the plot to the right
ax.legend(loc='upper left', bbox_to_anchor=(1.25, 0.75))

# Adjust the layout to make room for the legend
plt.subplots_adjust(right=0.5)



# ---------------------------------------------------------------------------- #
# Plot inlet states on p-s diagram
# ---------------------------------------------------------------------------- #
fig, ax = plt.subplots(figsize=(6.8, 4.8))
ax.set_xlabel("Entropy [kJ/kg/K]")
ax.set_ylabel("Pressure [bar]")
range_x = np.linspace(1000, 2000, 40)
range_y = np.linspace(59, 120, 40)*1e5
range_z = np.asarray([30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]) + 273.15
ax.set_xlim(np.asarray([range_x[0], range_x[-1]]) / 1e3)
ax.set_ylim(np.asarray([range_y[0], range_y[-1]]) / 1e5)
ax = sco2.plot_ps_diagram(fluid, ax, range_x, range_y, range_z)

# Loop over the unique values of 'tag'
for i, tag in enumerate(sorted(df["tag"].unique())):
    # Extract the subset of data for each tag.
    subset = df[df["tag"] == tag]

    # Define the label
    label = f"Dataset {tag}"
    
    # Plot the subset on the first subplot
    ax.plot(
        subset["s_in"],
        subset["p_in"],
        marker=markers[i % len(markers)],
        markersize=4.5,
        markerfacecolor="w",
        linewidth=0.75,
        linestyle="none",
        markeredgewidth=1.00,
        color=colors[i % len(colors)],
        label=label,
    )

# Place the legend outside of the figure
ax.legend(loc="upper right", fontsize=8)
fig.tight_layout(pad=1)

# Scale to nice units
bpy.scale_graphics_x(fig, +1e-3, mode="multiply")
bpy.scale_graphics_y(fig, +1e-5, mode="multiply")

# Save figure
if SAVE_FIGS:
    base_filename = os.path.join(OUTPUT_DIR, f"{CASE}_ps_diagram")
    bpy.savefig_in_formats(fig, base_filename, formats=[".png", ".svg"])


# ---------------------------------------------------------------------------- #
# Plot inlet states on T-s diagram
# ---------------------------------------------------------------------------- #
# Create figure
fig, ax = plt.subplots(figsize=(6.8, 4.8))
ax.set_xlabel("Entropy [kJ/kg/K]")
ax.set_ylabel("Temperature [$^\circ$C]")
range_x = np.linspace(1000, 2000, 40)
range_y = np.linspace(5, 55, 40) + 273.15
range_z = np.asarray([70, 75, 80, 85, 90, 95, 100]) * 1e5
ax.set_xlim(np.asarray([range_x[0], range_x[-1]]) / 1e3)
ax.set_ylim(np.asarray([range_y[0], range_y[-1]]) - 273.15)

# Plot diagram
ax = sco2.plot_Ts_diagram(fluid, ax, range_x, range_y, range_z)

# Loop over the unique values of 'tag'
for i, tag in enumerate(sorted(df["tag"].unique())):
    # Extract the subset of data for each tag.
    subset = df[df["tag"] == tag]

    # Define the label
    label = f"Dataset {tag}"
    
    # Plot the subset on the first subplot
    ax.plot(
        subset["s_in"],
        subset["T_in"],
        marker=markers[i % len(markers)],
        markersize=4.5,
        markerfacecolor="w",
        linewidth=0.75,
        linestyle="none",
        markeredgewidth=1.00,
        color=colors[i % len(colors)],
        label=label,
    )

# Place the legend outside of the figure
ax.legend(loc="lower right", fontsize=8)
fig.tight_layout(pad=1)

# Scale to nice units
bpy.scale_graphics_x(fig, +1e-3, mode="multiply")
bpy.scale_graphics_y(fig, -273.15, mode="add")

# Save figure
if SAVE_FIGS:
    base_filename = os.path.join(OUTPUT_DIR, f"{CASE}_Ts_diagram")
    bpy.savefig_in_formats(fig, base_filename, formats=[".png", ".svg"])


# Show the figures
plt.show()




