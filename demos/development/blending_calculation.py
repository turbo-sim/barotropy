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
fluid = bpy.Fluid(name="water", backend="HEOS")

# Create figure
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18.0, 5.0), gridspec_kw={"wspace": 0.25})
ax1.set_xlabel("Entropy (J/kg/K)")
ax1.set_ylabel("Temperature (K)")
ax2.set_xlabel("Pressure (Pa)")
ax2.set_ylabel("Density (kg/m$^3$)")
ax3.set_xlabel("Pressure (Pa)")
ax3.set_ylabel("Vapor quality (-)")
prop_x1, prop_y1 = "s","T"
prop_x2, prop_y2 = "p", "rho"
prop_x3, prop_y3 = "p", "Q"


# Plot phase diagram
fluid.plot_phase_diagram(
    prop_x1,
    prop_y1,
    axes=ax1,
    plot_critical_point=True,
    plot_quality_isolines=True,
    plot_pseudocritical_line=False,
    # plot_spinodal_line=True,
)

# Define inlet state
p_in, dT_subcooling = 60e5, 2
# p_in, dT_subcooling = 60e5, -2
props_in = bpy.compute_properties_coolprop(fluid._AS, bpy.PQ_INPUTS, p_in, 0)
props_in = bpy.compute_properties_coolprop(
    fluid.abstract_state,
    bpy.PT_INPUTS,
    p_in,
    props_in["T"] - dT_subcooling,
)

# Define phase change process
if props_in["s"] > fluid.critical_point.s:
    phase_change = "condensation"
    Q_onset = 0.99
    Q_width = 0.01
else:
    phase_change = "evaporation"
    Q_onset = 0.01
    Q_width = 0.01

# Define initial guess
rhoT_guess_metastable = [props_in["rho"], props_in["T"]]
rhoT_guess_equilibrium = [props_in["rho"], props_in["T"]]

# Loop over different exit pressures
p = np.linspace(p_in, p_in / 1.3, 50)
for i, p_out in enumerate(p):

    # Compute properties
    props_blend, props_eq, props_meta = bpy.compute_properties(
        fluid._AS,
        prop_1="p",
        prop_1_value=p_out,
        prop_2="s",
        prop_2_value=props_in["s"],
        rhoT_guess_equilibrium=rhoT_guess_equilibrium,
        rhoT_guess_metastable=rhoT_guess_metastable,
        calculation_type="blending",
        phase_change=phase_change,
        blending_variable="Q",
        blending_onset=Q_onset,
        blending_width=Q_width,
        generalize_quality=True,
        supersaturation=True,
        print_convergence=False,
        solver_algorithm="lm",
    )

    # Update initial guess
    rhoT_guess_metastable = [props_meta["rho"], props_meta["T"]]
    rhoT_guess_equilibrium = [props_eq["rho"], props_eq["T"]]

    # Make plots
    ax1.plot(props_eq[prop_x1], props_eq[prop_y1], marker="o", color=colors[0])
    ax1.plot(props_meta[prop_x1], props_meta[prop_y1], marker="o", color=colors[1])
    ax1.plot(props_blend[prop_x1], props_blend[prop_y1], marker="o", color=colors[3])

    ax2.plot(props_eq[prop_x2], props_eq[prop_y2], marker="o", color=colors[0])
    ax2.plot(props_meta[prop_x2], props_meta[prop_y2], marker="o", color=colors[1])
    ax2.plot(props_blend[prop_x2], props_blend[prop_y2], marker="o", color=colors[3])
    
    ax3.plot(props_eq[prop_x3], props_eq[prop_y3], marker="o", color=colors[0])
    ax3.plot(props_meta[prop_x3], props_meta[prop_y3], marker="o", color=colors[1])
    ax3.plot(props_blend[prop_x3], props_blend[prop_y3], marker="o", color=colors[3])


# Show figure
plt.show()
