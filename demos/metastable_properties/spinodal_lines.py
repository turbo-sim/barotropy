import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import barotropy as bpy

# Create the folder to save figures
bpy.set_plot_options(grid=False)
fig_dir = "output"
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

# Create fluid
fluid = bpy.Fluid(
    name="CO2",
    backend="HEOS",
    exceptions=True,
    generalize_quality=False,
    compute_subcooling=False,
    compute_superheating=False,
)

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
ax1.set_xlabel("Density (kg/m$^3$)")
ax1.set_ylabel("Pressure (Pa)")
ax2.set_xlabel("Density (kg/m$^3$)")
ax2.set_ylabel("Isothermal bulk modulus (Pa)")
ax1.set_ylim([fluid.triple_point_liquid.p, 2 * fluid.critical_point.p])
ax2.set_ylim([-2*fluid.critical_point.p, 2 * fluid.critical_point.p])
prop_x = "rhomass"
prop_y = "p"

# Create entropy range
rho_1 = fluid.triple_point_vapor.rho
rho_2 = fluid.triple_point_liquid.rho
delta_rho = rho_2 - rho_1
rho_array = np.logspace(np.log10(rho_1), np.log10(rho_2), 300)

# Compute equilibrium and metastable states
T_array = fluid.critical_point.T - np.asarray([5, 10, 15, 20, 25, 30])
states_eq = bpy.compute_property_grid(fluid, bpy.DmassT_INPUTS, rho_array, T_array)
states_meta = bpy.compute_property_grid_rhoT(fluid, rho_array, T_array)

# Plot properties
colormap = cm.magma(np.linspace(0.7, 0.1, len(T_array)))
for i, T in enumerate(T_array):
    # ax1.plot(
    #     states_eq[prop_x][i, :],
    #     states_eq[prop_y][i, :],
    #     color=colormap[i],
    #     linestyle="--"
    # )
    ax1.plot(
        states_meta[prop_x][i, :],
        states_meta[prop_y][i, :],
        color=colormap[i],
        label=f"$\Delta T_{{crit}}={fluid.critical_point.T-T:0.0f}$ K",
    )
    ax2.plot(
        states_meta[prop_x][i, :],
        states_meta["isothermal_bulk_modulus"][i, :],
        color=colormap[i],
        label=f"$\Delta T_{{crit}}={fluid.critical_point.T-T:0.0f}$ K",
    )


for i, T in enumerate(T_array):
    props = bpy.compute_spinodal_point(T, fluid._AS, branch="liquid", N_trial=100)
    ax1.plot(props[prop_x], props[prop_y], 
             marker="o", 
             markerfacecolor="w",
             color=colormap[i],)
    ax2.plot(props[prop_x], props["isothermal_bulk_modulus"], 
             marker="o", 
             markerfacecolor="w",
             color=colormap[i],)


    props = bpy.compute_spinodal_point(T, fluid._AS, branch="vapor", N_trial=100)
    ax1.plot(props[prop_x], props[prop_y],
             marker="o",
             markerfacecolor="w",
             color=colormap[i],)
    ax2.plot(props[prop_x], props["isothermal_bulk_modulus"], 
             marker="o",
             markerfacecolor="w",
             color=colormap[i],)


# Plot phase diagram
fluid.plot_phase_diagram(
    prop_x,
    prop_y,
    axes=ax1,
    plot_critical_point=True,
    plot_quality_isolines=False,
    plot_pseudocritical_line=False,
)
ax1.legend(loc="upper left", fontsize=10)
ax2.legend(loc="upper right", fontsize=10)
fig.tight_layout(pad=2)
bpy.savefig_in_formats(fig, os.path.join(fig_dir, f"spinodal_points_{fluid.name}"))

# Show figures
plt.show()



