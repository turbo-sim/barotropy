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


# State calculation using coolprop
fluid = bpy.Fluid(name="CO2", backend="HEOS")


# Create figure
fig, ax = plt.subplots(figsize=(6.0, 5.0))
ax.set_xlabel("Entropy (J/kg/K)")
ax.set_ylabel("Temperature (K)")
prop_x = "s"
prop_y = "T"


fig1, ax1 = plt.subplots(figsize=(6.0, 5.0))
ax1.set_xlabel("Pressure")
ax1.set_ylabel("Density")
prop_x2 = "p"
prop_y2 = "rho"

# Plot phase diagram
fluid.plot_phase_diagram(
    prop_x,
    prop_y,
    axes=ax,
    plot_critical_point=True,
    plot_quality_isolines=False,
    plot_pseudocritical_line=False,
    # plot_spinodal_line=True
)


# Compute subcooled inlet state
p_in, dT_subcooling = 60e5, 5



props_in = bpy.compute_properties_coolprop(fluid._AS, bpy.PQ_INPUTS, p_in, 0)
props_in = bpy.compute_properties_coolprop(
    fluid._AS, bpy.PT_INPUTS, p_in, props_in["T"] - dT_subcooling
)


s_in = props_in["s"]

p = np.linspace(p_in, p_in/2, 50)
for i, p_out in enumerate(p):


    props_out = bpy.compute_properties_coolprop(fluid._AS, bpy.PSmass_INPUTS, p_out, s_in)


    print(i, p_out)



    rho_guess1=props_in["rho"]
    T_guess1=props_in["T"]

    rhoT_guess_metastable = [rho_guess1, T_guess1]
    
    
    rho_guess=props_out["rho"]
    T_guess=props_out["T"]
    rhoT_guess_equilibrium = [rho_guess, T_guess]
    
    
    props_out = bpy.compute_properties(
        fluid._AS,
        prop_1="p",
        prop_1_value=p_out,
        prop_2="s",
        prop_2_value=s_in,
        rhoT_guess_equilibrium=rhoT_guess_equilibrium,
        calculation_type="equilibrium",
        print_convergence=False,
    )


    props_meta = bpy.compute_properties(
        fluid._AS,
        prop_1="p",
        prop_1_value=p_out,
        prop_2="s",
        prop_2_value=s_in,
        rhoT_guess_metastable=rhoT_guess_metastable,
        calculation_type="metastable",
        print_convergence=True,
    )


    props_blend = bpy.compute_properties(
        fluid._AS,
        prop_1="p",
        prop_1_value=p_out,
        prop_2="s",
        prop_2_value=s_in,
        rhoT_guess_equilibrium=rhoT_guess_equilibrium,
        rhoT_guess_metastable=rhoT_guess_metastable,
        calculation_type="blending",
        blending_variable="Q",
        blending_onset=0.00,
        blending_width=0.15,
        initial_phase="liquid",
        generalize_quality=True,
        supersaturation=True,
        print_convergence=False,
    )


    ax.plot(props_in[prop_x], props_in[prop_y], marker="o", label="in")
    ax.plot(props_out[prop_x], props_out[prop_y], marker="o", label="eq")
    # ax.plot(props_meta[prop_x], props_meta[prop_y], marker="o", label="meta")
    ax.plot(props_blend[prop_x], props_blend[prop_y], marker="+", label="blend")


    # ax.legend(loc="upper left")

    print(props_out["Q"])

    ax1.plot(props_in[prop_x2], props_in[prop_y2], marker="o", color="red")
    ax1.plot(props_out[prop_x2], props_out[prop_y2], marker="o", color="blue")
    ax1.plot(props_meta[prop_x2], props_meta[prop_y2], marker="o", color="green")
    ax1.plot(props_blend[prop_x2], props_blend[prop_y2], marker="+", color="black")




plt.show()
