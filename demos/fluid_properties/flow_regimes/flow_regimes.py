import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

import barotropy as bpy


# Create the folder to save figures
bpy.set_plot_options(grid=False)
colors = bpy.COLORS_MATLAB
fig_dir = "output"
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

# Define fluid
fluid = bpy.Fluid("CO2", backend="HEOS")
N_points = 100

# -------------------------------------------------------------------- #
# Compute properties along isentrope
# -------------------------------------------------------------------- #

# Compute subcooled inlet state
p_in, dT_subcooling = 0.75*fluid.critical_point.p, 5
p_out = p_in/2
Q_onset = 0.02
dQ_transition = 0.03
state_in = fluid.get_state(bpy.PQ_INPUTS, p_in, 0)
state_in = fluid.get_state(bpy.PT_INPUTS, p_in, state_in.T - dT_subcooling)

# Define initial guess
rhoT_guess_equilibrium = [state_in.rho, state_in.T]
rhoT_guess_metastable = [state_in.rho, state_in.T]

# Compute equilibrium, metastable and blending states
states_equilibrium, states_metastable, states_blended = [], [], []
p_array = np.linspace(p_in, p_out, N_points)
for p in p_array:
    state_blended, state_equilibrium, state_metastable = fluid.get_state_blending(
        prop_1="p",
        prop_1_value=p,
        prop_2="s",
        prop_2_value=state_in.s,
        rhoT_guess_equilibrium=rhoT_guess_equilibrium,
        rhoT_guess_metastable=rhoT_guess_metastable,
        phase_change="flashing",
        blending_variable="Q",
        blending_onset=Q_onset,
        blending_width=dQ_transition,
        print_convergence=False,
        supersaturation=True,
    )

    states_metastable.append(state_metastable)
    states_blended.append(state_blended)
    states_equilibrium.append(state_equilibrium)

    # Update initial guess
    rhoT_guess_equilibrium = [state_equilibrium.rho, state_equilibrium.T]
    rhoT_guess_metastable = [state_metastable.rho, state_metastable.T]

# Save last state
state_out_eq = states_equilibrium[-1]
state_out_meta = states_metastable[-1]
state_out_blend = states_blended[-1]

# Convert lists to dicts of arrays
states_equilibrium = bpy.states_to_dict(states_equilibrium)
states_metastable = bpy.states_to_dict(states_metastable)
states_blended = bpy.states_to_dict(states_blended)

# -------------------------------------------------------------------- #
# Plot phase diagram and void fraction contours
# -------------------------------------------------------------------- #

# Create a figure with two subplots side by side
fig, ax = plt.subplots(figsize=(6, 4.2))
prop_x = "s"
prop_y = "p"
prop_z = "void_fraction"
ax.set_xlabel(r"$s/s_\text{crit}$ - Reduced entropy")
ax.set_ylabel(r"$p/p_\text{crit}$ - Reduced pressure")
# ax.set_xlim([0.9 * fluid.triple_point_liquid[prop_x], 1.0 * fluid.triple_point_vapor[prop_x]])
# ax.set_ylim([1.2 * fluid.triple_point_liquid[prop_y], 1.1 * fluid.critical_point[prop_y]])
ax.set_xlim([0.25, 1.45])
ax.set_ylim([0.15, 1.25])

# Plots contour of void fraction
prop_dict = bpy.compute_quality_grid(fluid, num_points=100, quality_levels=np.linspace(0.0, 1.0, 100))
range_z = [0, 0.20, 0.80, 1]
# range_z = np.linspace(np.min(prop_dict[prop_z]), np.max(prop_dict[prop_z]), 10)
blues = plt.cm.Blues(np.linspace(0.4, 0.9, 256)[::-1])
colormap = mcolors.LinearSegmentedColormap.from_list("truncated_blues", blues)
contour = ax.contourf(
    prop_dict[prop_x],
    prop_dict[prop_y],
    prop_dict[prop_z],
    range_z,
    cmap=colormap,
)

# Add contour lines with black edges
contour_lines = ax.contour(
    prop_dict[prop_x],
    prop_dict[prop_y],
    prop_dict[prop_z],
    range_z,
    colors='black',
    linewidths=0.5,
)

# Add phase diagram
fluid.plot_phase_diagram(
    x_prop=prop_x,
    y_prop=prop_y,
    plot_saturation_line=True,
    # plot_spinodal_line=True,
    # plot_quality_isolines=True,
    N=N_points,
    axes=ax,
)

# Gather the indices for the first and last points
indices = [0, -1]
indices.append(np.argmin(np.abs(states_equilibrium["Q"])))
indices.append(np.argmin(np.abs(states_equilibrium["Q"] - Q_onset)))

# Plot all the selected points in a single command
ax.plot(
    states_blended[prop_x][indices],
    states_blended[prop_y][indices],
    color="black",
    linestyle="-",
    linewidth=1.00,
    marker="o",
    markersize=4,
)

# Scale plot
bpy.scale_graphics_x(fig, scale=1/fluid.critical_point[prop_x])
bpy.scale_graphics_y(fig, scale=1/fluid.critical_point[prop_y])



# -------------------------------------------------------------------- #
# Zoom-in axes
# -------------------------------------------------------------------- #

# Create the inset axis with more control over its position
# bbox_to_anchor=(left, bottom, width, height)
axins = inset_axes(ax, width="100%", height="100%", loc='upper left', 
                   bbox_to_anchor=(0.1, 0.6, 0.25, 0.35), bbox_transform=ax.transAxes)


# Reduce edge width of the inset axes
axins.tick_params(axis='both', which='both', width=0.75, length=3)  # Adjust width and length as needed
for spine in axins.spines.values():
    spine.set_linewidth(0.75)  # Adjust the edge width as needed

# Specify the part of the graph you want to zoom in on
x1, x2, y1, y2 = 0.73, 0.82, 0.56, 0.72  # adjust these values to focus on the desired area
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.contourf(
    prop_dict[prop_x]/fluid.critical_point[prop_x],
    prop_dict[prop_y]/fluid.critical_point[prop_y],
    prop_dict[prop_z],
    range_z,
    cmap=colormap,
)

# Add contour lines with black edges
axins.contour(
    prop_dict[prop_x]/fluid.critical_point[prop_x],
    prop_dict[prop_y]/fluid.critical_point[prop_y],
    prop_dict[prop_z],
    range_z,
    colors='black',
    linewidths=0.5,
)

# Plot all the selected points in a single command
axins.plot(
    states_blended[prop_x][indices]/fluid.critical_point[prop_x],
    states_blended[prop_y][indices]/fluid.critical_point[prop_y],
    color="black",
    linestyle="-",
    linewidth=1.00,
    marker="o",
    markersize=4.0,
)

# Create square around zoom axes
axins.tick_params(axis='both', which='major', labelsize=10)  # Adjust labelsize as needed
mark_inset(ax, axins, loc1=1, loc2=1, fc="none", ec="0", lw=0.75)


# Adjust pad
fig.tight_layout(pad=1)

# Save figure
bpy.savefig_in_formats(
    fig, os.path.join(fig_dir, f"flow_regimes_{fluid.name}")
)

# Show figures
plt.show()

