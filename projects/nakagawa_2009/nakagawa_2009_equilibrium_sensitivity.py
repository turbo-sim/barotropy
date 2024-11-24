import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import barotropy as bpy

bpy.print_package_info()
bpy.set_plot_options(grid=False)


# Create folder to save results
DIR_OUT = "./output/equilibrium_sensitivity"
os.makedirs(DIR_OUT, exist_ok=True)

# Barotropic figures settings
SHOW_FIG = False
SAVE_FIG = True
prop_x = "s"
prop_y = "p"

# Read case summary
DATAFILE = "./nakagawa_2009_cases.xlsx"
case_data = pd.read_excel(DATAFILE)

# Select cases of interest
CASES = ["nozzle_2"]
case_data = case_data[case_data["nozzle"].isin(CASES)]
case_data = case_data[case_data["tag"] == "base"]
# case_data = case_data[np.isclose(case_data["p0_in"], 61e5)]

# Create fluid
fluid_name = "CO2"
fluid = bpy.Fluid(name=fluid_name, backend="HEOS")

# Thermodynamic diagram
fig1, ax1 = plt.subplots(figsize=(6, 5))
ax1.set_xlabel(bpy.LABEL_MAPPING.get(prop_x, prop_x))
ax1.set_ylabel(bpy.LABEL_MAPPING.get(prop_y, prop_y))
fluid.plot_phase_diagram(
    prop_x, prop_y, axes=ax1, plot_quality_isolines=True, plot_spinodal_line=True
)
ax1.set_ylim([fluid.triple_point_liquid.p, 1.5 * fluid.critical_point.p])

# Efficiency sensitivity analysis
q_width = 0.04
q_onset_array = [0.00, 0.04, 0.08, 0.12]
# q_width = 0.05
# q_onset_array = [0.0, 0.05, 0.1, 0.15]
fig2, ax2 = plt.subplots(figsize=(6, 5))
ax2.set_xlabel(r"Pressure (Pa)")
ax2.set_ylabel(r"Density (kg/m$^3$)")
linestyles = ["-", "--", ":", "-."] * 3

# Loop over all cases
for i, (index, row) in enumerate(case_data.iterrows()):

    # Import experimental data and convert to SI units
    print(f"Processing case {index}")
    T_in = row["T0_in"] + 273.15
    p_in = row["p0_in"]
    # p_out = 0.50 * fluid.critical_point.p
    p_out = 1.05 * fluid.triple_point_liquid.p
    # polytropic_efficiency = row["polytropic_efficiency"]
    polytropic_efficiency = 0.2

    model_eq = bpy.BarotropicModel(
        fluid_name=fluid_name,
        T_in=T_in,
        p_in=p_in,
        p_out=p_out,
        efficiency=polytropic_efficiency,
        calculation_type="equilibrium",
        HEOS_solver="hybr",
        ODE_solver="LSODA",
        ODE_tolerance=1e-9,
        polynomial_degree=4,
        polynomial_format="horner",
        output_dir=DIR_OUT,
    )
    model_eq.solve()

    
    # Loop over efficiencies
    for j, q_onset in enumerate(q_onset_array):
        model = bpy.BarotropicModel(
            fluid_name=fluid_name,
            T_in=T_in,
            p_in=p_in,
            p_out=p_out,
            efficiency=polytropic_efficiency,
            calculation_type="blending",
            blending_onset=q_onset,
            blending_width=q_width,
            HEOS_solver="hybr",
            ODE_solver="LSODA",
            ODE_tolerance=1e-9,
            polynomial_degree=[4, 4, 8],
            polynomial_format="horner",
            output_dir=DIR_OUT,
        )
        model.solve()
        model.fit_polynomials()
        if j == 0:
            ax1.plot(
                model.states[prop_x],
                model.states[prop_y],
                label=rf"$p_\text{{in}}={{{p_in/1e5:0.1f}}} \, \mathrm{{bar}}$, $T_\text{{in}}={{{T_in-273.15:0.1f}}} \, ^\circ\mathrm{{C}}$",
                color=bpy.COLORS_MATLAB[i],
            )
            ax2.plot(
                model_eq.states["p"],
                model_eq.states["density"],
                color="black",
                linestyle="-",
            )

        ax2.plot(
            model.states["p"],
            model.states["density"],
            label=rf"$Q_\text{{onset}}={{{q_onset*100:0.0f}}} \, \%$",
            color=bpy.COLORS_MATLAB[i],
            linestyle=linestyles[j],
        )

# Save diagram
ax1.legend(loc="upper left", fontsize=11)
fig1.tight_layout(pad=1)
bpy.savefig_in_formats(fig1, os.path.join(DIR_OUT, "thermodynamic_diagram"))

# Save density plots
plt.legend(
    loc="lower right",
    fontsize=11,
    ncols=len(case_data),
    handletextpad=0.4,  # Space between legend marker and text
    columnspacing=0.5,  # Space between columns
    borderaxespad=0.2,  # Padding between legend and plot area
)
fig2.tight_layout(pad=1)
bpy.savefig_in_formats(fig2, os.path.join(DIR_OUT, f"equilibrium_sensitivity"))

# Show figure
plt.show()
