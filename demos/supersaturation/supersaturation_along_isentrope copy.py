import os
import numpy as np
import matplotlib.pyplot as plt

import barotropy as bpy


# Create the folder to save figures
bpy.set_plot_options(grid=False)
colors = bpy.COLORS_MATLAB
fig_dir = "output"
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

# Define fluid
fluid = bpy.Fluid("CO2", backend="HEOS")


# -------------------------------------------------------------------- #
# Compute degree of supersaturation along liquid isentrope
# -------------------------------------------------------------------- #

# Compute subcooled inlet state
p_in, dT_subcooling = 60e5, 5
p_out = 30e5
state_in = fluid.get_state(bpy.PQ_INPUTS, p_in, 0)
state_in = fluid.get_state(bpy.PT_INPUTS, p_in, state_in.T - dT_subcooling)

# Compute equilibrium exit states
p_array = np.linspace(p_in, p_out, 100)
states_eq = []
for p in p_array:
    states_eq.append(fluid.get_state(bpy.PSmass_INPUTS, p, state_in.s))
state_out_eq = states_eq[-1]
states_eq = bpy.states_to_dict(states_eq)

# Compute metastable exit states (use inlet state as initial guess)
states_meta = []
for p in p_array:
    states_meta.append(
        fluid.get_state_metastable(
            "p",
            p,
            "s",
            state_in.s,
            rho_guess=state_in.rho,
            T_guess=state_in.T,
            print_convergence=False,
            supersaturation=True,
        )
    )
state_out_meta = states_meta[-1]
states_meta = bpy.states_to_dict(states_meta)

# Print properties at the inlet and outlet
msgs = [
    f"{'State':>10} {'Temperature (K)':>20} {'Pressure (Pa)':>20} {'Density (kg/m3)':>20} {'Enthalpy (J/kg)':>20} {'Entropy (J/kg/K)':>20} {'Supersaturation (K)':>20} {'Supersaturation ratio':>20}",
    f"{'1':>10} {state_in.T:20.4f} {state_in.p:20.1f} {state_in.rho:20.1f} {state_in.h:20.1f} {state_in.s:20.1f} {np.nan:>20} {np.nan:>20}",
    f"{'2':>10} {state_out_eq.T:20.4f} {state_out_eq.p:20.1f} {state_out_eq.rho:20.1f} {state_out_eq.h:20.1f} {state_out_eq.s:20.1f} {np.nan:>20} {np.nan:>20}",
    f"{'2m':>10} {state_out_meta.T:20.4f} {state_out_meta.p:20.1f} {state_out_meta.rho:20.1f} {state_out_meta.h:20.1f} {state_out_meta.s:20.1f} {state_out_meta.supersaturation_degree:>20.03f} {state_out_meta.supersaturation_ratio:>20.03f}",
    "",
]
for msg in msgs:
    print(msg)

# Create a figure with two subplots side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"wspace": 0.25})

# Plot phase diagram and thermodynamic points on the first subplot
x1 = "s"
y1 = "T"
ax1.set_xlabel("Entropy (J/kg/K))")
ax1.set_ylabel("Temperature (K)")
ax1.set_xlim([1.5 * fluid.triple_point_liquid[x1], 0.9 * fluid.triple_point_vapor[x1]])
ax1.set_ylim([1.2 * fluid.triple_point_liquid[y1], 1.1 * fluid.critical_point[y1]])
fluid.plot_phase_diagram(
    x_prop=x1,
    y_prop=y1,
    plot_saturation_line=True,
    plot_spinodal_line=True,
    plot_quality_isolines=True,
    N=50,
    axes=ax1,
)

# Plot equilibrium and metastable states
ax1.plot(
    states_eq[x1][[0, -1]]*1.005,
    states_eq[y1][[0, -1]],
    color=colors[0],
    linestyle="-",
    marker="o",
    label="Equilibrium expansion",
)
ax1.plot(
    states_meta[x1][[0, -1]],
    states_meta[y1][[0, -1]],
    color=colors[1],
    linestyle="-",
    marker="o",
    label="Metastable expansion",
)
ax1.legend(loc="upper right")

# Define x and y variables for the second subplot
x2 = "p"
y2 = "rho"
y2_bis = "supersaturation_degree"

# Plot the equilibrium and metastable states on the primary y-axis of the second subplot
ax2.set_xlabel("Pressure (Pa)")
ax2.set_ylabel("Density (kg/m$^3$)")
ax2.set_ylim([0, 1.2*states_eq[y2][0]])
line1, = ax2.plot(states_eq[x2], states_eq[y2], color=colors[0], linestyle="-", label="Equilibrium density")
line2, = ax2.plot(states_meta[x2], states_meta[y2], color=colors[1], linestyle="-", label="Metastable density")

# Create a secondary y-axis for the supersaturation degree
ax2_secondary = ax2.twinx()
ax2_secondary.set_ylabel("Supersaturation degree (K)")
ax2_secondary.set_ylim(1.2*np.asarray([-1, 1])*np.max(states_meta[y2_bis]))
line3, = ax2_secondary.plot(states_meta[x2], states_meta[y2_bis], color=colors[1], linestyle="--", label="Supersaturation degree")

# Combine legends from both y-axes
lines = [line1, line2, line3]
labels = [line.get_label() for line in lines]
ax2.legend(lines, labels, loc="lower right")

bpy.savefig_in_formats(fig, os.path.join(fig_dir, f"supersaturation_along_liquid_isentrope_{fluid.name}"))



# -------------------------------------------------------------------- #
# Compute degree of supersaturation along vapor isentrope
# -------------------------------------------------------------------- #

# Compute subcooled inlet state
p_in, dT_superheating = 70e5, 5
p_out = 50e5
state_in = fluid.get_state(bpy.PQ_INPUTS, p_in, 0)
state_in = fluid.get_state(bpy.PT_INPUTS, p_in, state_in.T + dT_superheating)

# Compute equilibrium exit states
p_array = np.linspace(p_in, p_out, 100)
states_eq = []
for p in p_array:
    states_eq.append(fluid.get_state(bpy.PSmass_INPUTS, p, state_in.s))
state_out_eq = states_eq[-1]
states_eq = bpy.states_to_dict(states_eq)

# Compute metastable exit states (use inlet state as initial guess)
states_meta = []
for p in p_array:
    states_meta.append(
        fluid.get_state_metastable(
            "p",
            p,
            "s",
            state_in.s,
            rho_guess=state_in.rho,
            T_guess=state_in.T,
            print_convergence=False,
            supersaturation=True,
        )
    )
state_out_meta = states_meta[-1]
states_meta = bpy.states_to_dict(states_meta)

# Print properties at the inlet and outlet
msgs = [
    f"{'State':>10} {'Temperature (K)':>20} {'Pressure (Pa)':>20} {'Density (kg/m3)':>20} {'Enthalpy (J/kg)':>20} {'Entropy (J/kg/K)':>20} {'Supersaturation (K)':>20} {'Supersaturation ratio':>20}",
    f"{'1':>10} {state_in.T:20.4f} {state_in.p:20.1f} {state_in.rho:20.1f} {state_in.h:20.1f} {state_in.s:20.1f} {np.nan:>20} {np.nan:>20}",
    f"{'2':>10} {state_out_eq.T:20.4f} {state_out_eq.p:20.1f} {state_out_eq.rho:20.1f} {state_out_eq.h:20.1f} {state_out_eq.s:20.1f} {np.nan:>20} {np.nan:>20}",
    f"{'2m':>10} {state_out_meta.T:20.4f} {state_out_meta.p:20.1f} {state_out_meta.rho:20.1f} {state_out_meta.h:20.1f} {state_out_meta.s:20.1f} {state_out_meta.supersaturation_degree:>20.03f} {state_out_meta.supersaturation_ratio:>20.03f}",
    "",
]
for msg in msgs:
    print(msg)

# Create a figure with two subplots side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"wspace": 0.25})

# Plot phase diagram and thermodynamic points on the first subplot
x1 = "s"
y1 = "T"
ax1.set_xlabel("Entropy (J/kg/K))")
ax1.set_ylabel("Temperature (K)")
ax1.set_xlim([1.5 * fluid.triple_point_liquid[x1], 0.9 * fluid.triple_point_vapor[x1]])
ax1.set_ylim([1.2 * fluid.triple_point_liquid[y1], 1.1 * fluid.critical_point[y1]])
fluid.plot_phase_diagram(
    x_prop=x1,
    y_prop=y1,
    plot_saturation_line=True,
    plot_spinodal_line=True,
    plot_quality_isolines=True,
    N=50,
    axes=ax1,
)

# Plot equilibrium and metastable states
ax1.plot(
    states_meta[x1][[0, -1]],
    states_meta[y1][[0, -1]],
    color=colors[1],
    linestyle="-",
    marker="o",
    label="Metastable expansion",
)
ax1.plot(
    states_eq[x1][[0, -1]]*1.005,
    states_eq[y1][[0, -1]],
    color=colors[0],
    linestyle="-",
    marker="o",
    label="Equilibrium expansion",
)
ax1.legend(loc="upper right")

# Define x and y variables for the second subplot
x2 = "p"
y2 = "rho"
y2_bis = "supersaturation_degree"

# Plot the equilibrium and metastable states on the primary y-axis of the second subplot
ax2.set_xlabel("Pressure (Pa)")
ax2.set_ylabel("Density (kg/m$^3$)")
ax2.set_ylim([0, 1.2*states_eq[y2][0]])
line1, = ax2.plot(states_eq[x2], states_eq[y2], color=colors[0], linestyle="-", label="Equilibrium density")
line2, = ax2.plot(states_meta[x2], states_meta[y2], color=colors[1], linestyle="-", label="Metastable density")

# Create a secondary y-axis for the supersaturation degree
ax2_secondary = ax2.twinx()
ax2_secondary.set_ylabel("Supersaturation degree (K)")
ax2_secondary.set_ylim(1.4*np.asarray([-1, 1])*np.max(states_meta[y2_bis]))
line3, = ax2_secondary.plot(states_meta[x2], states_meta[y2_bis], color=colors[1], linestyle="--", label="Supersaturation degree")
# ax2.axhline(y=0, color="black")

# Combine legends from both y-axes
lines = [line1, line2, line3]
labels = [line.get_label() for line in lines]
ax2.legend(lines, labels, loc="lower right")

bpy.savefig_in_formats(fig, os.path.join(fig_dir, f"supersaturation_along_vapor_isentrope_{fluid.name}"))



# Show figures
plt.show()
