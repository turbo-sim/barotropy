import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import barotropy as bpy

# Create the folder to save figures
bpy.set_plot_options()
fig_dir = "output"
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

# Create fluid
fluid = bpy.Fluid(
    name="CO2",
    exceptions=True,
    generalize_quality=True,
    compute_subcooling=True,
    compute_superheating=True,
)



# Create figure
fig, ax = plt.subplots(figsize=(6.0, 5.0))
ax.set_xlabel("Entropy (J/kg/K)")
ax.set_ylabel("Temperature (K)")
prop_x = "s"
prop_y = "T"

# Plot phase diagram
fluid.plot_phase_diagram(
    prop_x,
    prop_y,
    axes=ax,
    plot_critical_point=True,
    plot_quality_isolines=False,
    plot_pseudocritical_line=False,
)


# Try a p-T function call
p = fluid.critical_point.p * 0.9
T = fluid.critical_point.T * 0.99
state = fluid.set_state(bpy.PT_INPUTS, p, T)
rho_guess = state.rho * 1.1
T_guess = state.T - 2

state_meta = fluid.set_state_metastable(
    prop_1="p",
    prop_1_value=p,
    prop_2="T",
    prop_2_value=T,
    rho_guess=rho_guess,
    T_guess=T_guess,
    print_convergence=True,
)

ax.plot(state[prop_x], state[prop_y], "o")
ax.plot(state_meta[prop_x], state_meta[prop_y], "+")



# # Try a trivial rho-T function call
# rho = 1.2*fluid.critical_point.rho
# T = 30 + fluid.critical_point.T
# state = fluid.set_state(bpy.DmassT_INPUTS, rho, T)
# state_meta = fluid.set_state_metastable(
#     prop_1="rhomass",
#     prop_1_value=rho,
#     prop_2="T",
#     prop_2_value=T,
#     rho_guess=rho,
#     T_guess=T,
#     print_convergence=True,
# )

# ax.plot(state[prop_x], state[prop_y], "o")
# ax.plot(state_meta[prop_x], state_meta[prop_y], "+")


plt.show()
