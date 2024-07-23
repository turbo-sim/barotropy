import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import barotropy as bpy


# Create the folder to save figures
bpy.set_plot_options()
fig_dir = "figures"
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

# Create fluid
fluid = bpy.Fluid(name="CO2", exceptions=True)


# --------------------------------------------------------------------------- #
# Plot iso-pressure lines
# --------------------------------------------------------------------------- #

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
ax1.set_xlabel("Entropy (J/kg/K)")
ax1.set_ylabel("Temperature (K)")
ax2.set_xlabel("Entropy (J/kg/K)")
ax2.set_ylabel("Vapor quality (-)")
# ax1.set_ylim([0.75 * fluid.critical_point.T, 1.25 * fluid.critical_point.T])
prop_x = "s"
prop_y = "T"

# Create entropy range
s1 = fluid.triple_point_liquid.s
s2 = fluid.triple_point_vapor.s
delta_s = s2 - s1
s_array = np.linspace(s1 + delta_s / 8, s2 + delta_s / 16, 100)

# Subcritical cases
p_array = np.asarray([0.5, 0.6, 0.7, 0.8, 0.9, 0.99]) * fluid.critical_point.p
states = bpy.compute_properties_meshgrid(fluid, bpy.PSmass_INPUTS, p_array, s_array, generalize_quality=True)
colormap = cm.magma(np.linspace(0.1, 0.7, len(p_array)))
for i in range(states[prop_x].shape[-1]):
    ax1.plot(
        states[prop_x][:, i],
        states[prop_y][:, i],
        color=colormap[i],
        label=f"$p/p_{{crit}}={p_array[i]/fluid.critical_point.p:0.2f}$ bar",
    )
    ax2.plot(
        states[prop_x][:, i],
        states["Q"][:, i],
        color=colormap[i],
        label=f"$p/p_{{crit}}={p_array[i]/fluid.critical_point.p:0.2f}$ bar",
    )

# Supercritical cases
p_array = np.asarray([1.01, 1.2, 1.4, 1.6, 1.8, 2.0]) * fluid.critical_point.p
states = bpy.compute_properties_meshgrid(fluid, bpy.PSmass_INPUTS, p_array, s_array, generalize_quality=True)
colormap = cm.magma(np.linspace(0.7, 0.1, len(p_array)))
for i in range(states[prop_x].shape[-1]):
    ax1.plot(
        states[prop_x][:, i],
        states[prop_y][:, i],
        color=colormap[i],
        linestyle="--",
        label=f"$p/p_{{crit}}={p_array[i]/fluid.critical_point.p:0.2f}$ bar",
    )
    ax2.plot(
        states[prop_x][:, i],
        states["Q"][:, i],
        color=colormap[i],
        linestyle="--",
        label=f"$p/p_{{crit}}={p_array[i]/fluid.critical_point.p:0.2f}$ bar",
    )

# Plot phase diagram
fluid.plot_phase_diagram(
    prop_x,
    prop_y,
    axes=ax1,
    plot_critical_point=True,
    plot_quality_isolines=True,
    plot_pseudocritical_line=True,
)
ax1.legend(loc="upper left", fontsize=10)
ax2.legend(loc="upper left", fontsize=10)
fig.tight_layout(pad=2)
bpy.savefig_in_formats(fig, os.path.join(fig_dir, "extended_quality_isobars"))

# p_array1 = np.asarray(np.linspace(0.5, 0.99, 100)) * fluid.critical_point.p
# p_array2 = np.asarray(np.linspace(1.01, 2.00, 100)) * fluid.critical_point.p
# p_array = np.concatenate([p_array1, p_array2])
# states = bpy.compute_properties_meshgrid(fluid, bpy.PSmass_INPUTS, p_array, s_array, generalize_quality=True)
# contour = ax1.contour(
#     states[prop_x],
#     states[prop_y],
#     states["Q"],
#     np.linspace(-1, 2, 31),
#     linewidths=0.5,
# )


# --------------------------------------------------------------------------- #
# Plot iso-temperature lines
# --------------------------------------------------------------------------- #

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
ax1.set_xlabel("Entropy (J/kg/K)")
ax1.set_ylabel("Temperature (K)")
ax2.set_xlabel("Entropy (J/kg/K)")
ax2.set_ylabel("Vapor quality (-)")
# ax1.set_ylim([0.5 * fluid.critical_point.T, 1.5 * fluid.critical_point.T])
prop_x = "s"
prop_y = "T"

# Create entropy range
s1 = fluid.triple_point_liquid.s
s2 = fluid.triple_point_vapor.s
delta_s = s2 - s1
s_array = np.linspace(s1 + delta_s / 8, s2 + delta_s / 16, 100)

# Subcritical cases
T_array = np.asarray([0.75, 0.8, 0.85, 0.9, 0.95, 0.99]) * fluid.critical_point.T
colormap = cm.magma(np.linspace(0.1, 0.7, len(T_array)))
states = bpy.compute_properties_meshgrid(fluid, bpy.SmassT_INPUTS, s_array, T_array, generalize_quality=True)
for i in range(states[prop_x].shape[0]):
    ax1.plot(
        states[prop_x][i, :],
        states[prop_y][i, :],
        color=colormap[i],
        label=f"$T/T_{{crit}}={T_array[i]/fluid.critical_point.p:0.2f}$ bar",
    )
    ax2.plot(
        states[prop_x][i, :],
        states["Q"][i, :],
        color=colormap[i],
        label=f"$T/T_{{crit}}={T_array[i]/fluid.critical_point.p:0.2f}$ bar",
    )

# Supercritical cases
T_array = np.asarray([1.01, 1.1, 1.2, 1.3, 1.4, 1.5]) * fluid.critical_point.T
colormap = cm.magma(np.linspace(0.7, 0.1, len(p_array)))
states = bpy.compute_properties_meshgrid(fluid, bpy.SmassT_INPUTS, s_array, T_array, generalize_quality=True)
for i in range(states[prop_x].shape[0]):
    ax1.plot(
        states[prop_x][i, :],
        states[prop_y][i, :],
        color=colormap[i],
        linestyle="--",
        label=f"$T/T_{{crit}}={T_array[i]/fluid.critical_point.p:0.2f}$ bar",
    )
    ax2.plot(
        states[prop_x][i, :],
        states["Q"][i, :],
        color=colormap[i],
        linestyle="--",
        label=f"$T/T_{{crit}}={T_array[i]/fluid.critical_point.p:0.2f}$ bar",
    )

# Plot phase diagram
fluid.plot_phase_diagram(
    prop_x,
    prop_y,
    axes=ax1,
    plot_critical_point=True,
    plot_quality_isolines=True,
    plot_pseudocritical_line=True,
)
ax1.legend(loc="upper left", fontsize=10)
ax2.legend(loc="upper left", fontsize=10)
fig.tight_layout(pad=2)
bpy.savefig_in_formats(fig, os.path.join(fig_dir, "extended_quality_isotherms"))


# --------------------------------------------------------------------------- #
# Plot iso-entropy lines
# --------------------------------------------------------------------------- #

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
ax1.set_xlabel("Entropy (J/kg/K)")
ax1.set_ylabel("Temperature (K)")
ax2.invert_xaxis()
ax2.set_xlabel("Temperature (K)")
ax2.set_ylabel("Vapor quality (-)")
# ax1.set_ylim([0.5 * fluid.critical_point.T, 1.5 * fluid.critical_point.T])
prop_x = "s"
prop_y = "T"

# Create temperature range
T_array = np.linspace(0.75, 1.5, 101) * fluid.critical_point.T

# Liquid-like cases
s_array = np.asarray([0.7, 0.8, 0.9, 0.99]) * fluid.critical_point.s
colormap = cm.magma(np.linspace(0.1, 0.7, len(s_array)))
states = bpy.compute_properties_meshgrid(fluid, bpy.SmassT_INPUTS, s_array, T_array, generalize_quality=True)
for i in range(states[prop_x].shape[1]):
    ax1.plot(
        states[prop_x][:, i],
        states[prop_y][:, i],
        color=colormap[i],
        label=f"$s/s_{{crit}}={s_array[i]/fluid.critical_point.s:0.2f}$ bar",
    )
    ax2.plot(
        states[prop_y][:, i],
        states["Q"][:, i],
        color=colormap[i],
        label=f"$s/s_{{crit}}={s_array[i]/fluid.critical_point.s:0.2f}$ bar",
    )

# Gas-like cases
s_array = np.asarray([1.01, 1.1, 1.2, 1.3]) * fluid.critical_point.s
colormap = cm.magma(np.linspace(0.7, 0.1, len(s_array)))
states = bpy.compute_properties_meshgrid(fluid, bpy.SmassT_INPUTS, s_array, T_array, generalize_quality=True)
for i in range(states[prop_x].shape[1]):
    print(states[prop_x][:, i])
    ax1.plot(
        states[prop_x][:, i],
        states[prop_y][:, i],
        color=colormap[i],
        linestyle="--",
        label=f"$s/s_{{crit}}={s_array[i]/fluid.critical_point.s:0.2f}$ bar",
    )
    ax2.plot(
        states[prop_y][:, i],
        states["Q"][:, i],
        color=colormap[i],
        linestyle="--",
        label=f"$s/s_{{crit}}={s_array[i]/fluid.critical_point.s:0.2f}$ bar",
    )

# Plot phase diagram
fluid.plot_phase_diagram(
    prop_x,
    prop_y,
    axes=ax1,
    plot_critical_point=True,
    plot_quality_isolines=True,
    plot_pseudocritical_line=True,
)
ax1.legend(loc="upper left", fontsize=10)
ax2.legend(loc="upper left", fontsize=10)
fig.tight_layout(pad=2)
bpy.savefig_in_formats(fig, os.path.join(fig_dir, "extended_quality_isentropes"))


# Show figures
plt.show()
