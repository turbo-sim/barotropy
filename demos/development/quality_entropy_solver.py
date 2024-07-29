import numpy as np
import matplotlib.pyplot as plt

import barotropy as bpy
import barotropy.fluids.core as props

# Observations:
# This inverse call is not perfect, it fails sometimes close to critical point
# Overall the convergence is quite good
# Note that there are some points with non-unique Q-s solutions close to the Q=0.5 line

# Create fluid
fluid_name = "CO2"
fluid = bpy.Fluid(name=fluid_name)

# Arrays of pressures (Pa) and qualities
pressures = [50e5, 60e5, 70e5]*4  # Example pressures
qualities = [0.00]*3 + [0.20]*3 + [0.80]*3 + [1.00]*3  # Example qualities

# Create figure for plotting
fig, ax = plt.subplots(figsize=(6.8, 4.8))
ax.set_xlabel("Entropy [kJ/kg/K]")
ax.set_ylabel("Temperature [$^\circ$C]")

# Plot phase diagram
ax = fluid.plot_phase_diagram("s", "T", plot_quality_isolines=True)

# Iterate over the arrays of pressures and qualities
for p, Q in zip(pressures, qualities):
    # Set state based on current pressure and quality
    state_1 = fluid.get_state(props.PQ_INPUTS, p, Q)
    state_2 = props.set_state_Qs(fluid, Q, state_1.s)
    
    # Plot state
    ax.plot(state_1.s, state_1.T, marker="o", color=bpy.COLORS_MATLAB[0])
    ax.plot(state_2.s, state_2.T, marker="+", color=bpy.COLORS_MATLAB[1])
    
    # Print data
    print(f"Pressure: {state_1.p/1e5:0.6e} bar | Entropy: {state_1.s/1e3:0.6e} kJ/kg/K | Quality: {state_1.Q:0.6e}")
    print(f"Pressure: {state_2.p/1e5:0.6e} bar | Entropy: {state_2.s/1e3:0.6e} kJ/kg/K | Quality: {state_2.Q:0.6e}")
    print()

# Show figure
plt.show()

