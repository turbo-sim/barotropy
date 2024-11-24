import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import barotropy as bpy

bpy.print_package_info()
bpy.set_plot_options(grid=False)


# Create folder to save results
DIR_OUT = "./output/efficiency_sensitivity"
os.makedirs(DIR_OUT, exist_ok=True)

# Barotropic figures settings
SHOW_FIG = False
SAVE_FIG = True
prop_x = "s"
prop_y = "T"

# Read case summary
DATAFILE = "./lettieri_2018_data.xlsx"
case_data = pd.read_excel(DATAFILE)

# Select cases of interest
# CASES = ["nozzle_2"]
case_data = case_data.iloc[::2, :]
# case_data = case_data[case_data["nozzle"].isin(CASES)]
# case_data = case_data[case_data["tag"] == "a"]
# case_data = case_data[np.isclose(case_data["p0_in"], 61e5)]

# Create fluid
fluid_name = "CO2"
fluid = bpy.Fluid(name=fluid_name, backend="HEOS")

# Thermodynamic diagram 
fig1, ax1 = plt.subplots(figsize=(6, 5))
ax1.set_xlabel(bpy.LABEL_MAPPING.get(prop_x, prop_x))
ax1.set_ylabel(bpy.LABEL_MAPPING.get(prop_y, prop_y))
fluid.plot_phase_diagram(prop_x, prop_y, axes=ax1, plot_quality_isolines=True, plot_spinodal_line=False)

# Efficiency sensitivity analysis
efficiency_array = [1.00, 0.90, 0.80]
fig2, ax2 = plt.subplots(figsize=(6, 5))
ax2.set_xlabel(r"Pressure (Pa)")
ax2.set_ylabel(r"Density (kg/m$^3$)")
linestyles = ['-', '--', ':', '-.'] * 3
colors = plt.cm.magma(np.linspace(0.25, 0.75, len(case_data)))

# Loop over all cases
for i, (index, row) in enumerate(case_data.iterrows()):

    # Import experimental data and convert to SI units
    print(f"Processing case {index}")
    T_in = row["T0_in"] + 273.15
    p_in = row["p0_in"] * 1e5
    p_out = 1.1*fluid.triple_point_vapor.p
    polytropic_efficiency = row["polytropic_efficiency"]
    calculation_type = row["calculation_type"]
    q_onset = row["q_onset"]
    q_transition = row["q_transition"]

    # Create barotropic model object
    model = bpy.BarotropicModel(
        fluid_name=fluid_name,
        T_in=T_in,
        p_in=p_in,
        p_out=p_out,
        efficiency=polytropic_efficiency,
        calculation_type=calculation_type,
        blending_onset=q_onset,
        blending_width=q_transition,
        HEOS_solver="hybr",
        ODE_solver="LSODA",
        ODE_tolerance=1e-9,
        polynomial_degree=4,
        polynomial_format="horner",
        output_dir=DIR_OUT,
    )
    

    # Evaluate barotropic model and export polynomial expressions
    model.solve()


    # Loop over efficiencies
    for j, efficiency in enumerate(efficiency_array):
        model_sens = bpy.BarotropicModel(
            fluid_name=fluid_name,
            T_in=T_in,
            p_in=p_in,
            p_out=p_out,
            efficiency=efficiency,
            calculation_type="equilibrium",
            # blending_onset=q_onset,
            # blending_width=q_transition,
            HEOS_solver="hybr",
            ODE_solver="LSODA",
            ODE_tolerance=1e-9,
            polynomial_degree=4,
            polynomial_format="horner",
            output_dir=DIR_OUT,
        )
        model_sens.solve()
        ax1.plot(
            model_sens.states[prop_x],
            model_sens.states[prop_y],
            label=rf"$p_\text{{in}}={{{p_in/1e5:0.1f}}} \, \mathrm{{bar}}$, $T_\text{{in}}={{{T_in-273.15:0.1f}}} \, ^\circ\mathrm{{C}}$" if j == 1 else None,
            color=colors[i],
            linestyle=linestyles[j],
        )
        ax2.plot(
            model_sens.states["p"],
            model_sens.states["density"],
            label=rf"$\eta_\text{{p}}={{{efficiency*100:0.0f}}} \, \%$" if i == 1 else None,
            color=colors[i],
            linestyle=linestyles[j],
        )

# Save diagram
ax1.legend(loc="upper left", fontsize=12)
fig1.tight_layout(pad=1)
bpy.savefig_in_formats(fig1, os.path.join(DIR_OUT, "thermodynamic_diagram"))

# Save density plots
plt.legend(
    loc="upper left",
    fontsize=12,
    ncols=1,
    handletextpad=0.4,  # Space between legend marker and text
    columnspacing=0.5,  # Space between columns
    borderaxespad=0.2,  # Padding between legend and plot area
)
fig2.tight_layout(pad=1)
bpy.savefig_in_formats(fig2, os.path.join(DIR_OUT, f"efficiency_sensitivity"))

# Show figure
plt.show()


