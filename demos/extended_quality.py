
import numpy as np
import matplotlib.pyplot as plt

import barotropy as bpy
import barotropy.fluids.core as props

fluid_name = "CO2"

fluid = bpy.Fluid(name=fluid_name)


p = 50e5
liq_sat = fluid.set_state(props.PQ_INPUTS, p, 0.0)
vap_sat = fluid.set_state(props.PQ_INPUTS, p, 1.0)
h_array = np.linspace(0.8*liq_sat.h, 1.2*vap_sat.h, 250)
states = []

p = 190e5
for h in h_array:
    states.append(fluid.set_state(props.HmassP_INPUTS, h, p))
properties = props.states_to_dict(states)

print(properties["Q"])

# Plot evolution of flow variables
bpy.set_plot_options()
figure, ax = plt.subplots(figsize=(6.0, 4.8))
# ax.set_xlabel("Entropy")
ax.set_ylabel("Quality")
ax.plot(
    properties["s"],
    properties["Q"],
    linewidth=1.00,
    marker="o",
    markersize=3.5,
    markeredgewidth=1.00,
    markerfacecolor="w",
)


plt.show()