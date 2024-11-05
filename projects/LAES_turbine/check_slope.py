import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import barotropy as bpy

bpy.print_package_info()
bpy.set_plot_options()


from matplotlib import cm

# Create folder to save results
# DIR_OUT = "./output"
DIR_OUT = "./barotropic_model_20240805"
if not os.path.isdir(DIR_OUT):
    os.makedirs(DIR_OUT)

# Load cases from Excel
case_table = pd.read_excel("./simulation_cases.xlsx")
case_indices = case_table["index"]
# case_indices = [1, 2]
case_indices = [2]
case_table = case_table[case_table["index"].isin(case_indices)]

# Define plot settings
var_x = "s"
var_y = "T"
save_figures = True
colors = cm.magma(np.linspace(0.2, 0.8, max(len(case_indices) - 1 , 1)))  # Reds, Blues, magma

# Loop over cases
for i, (idx, row) in enumerate(case_table.iterrows()):

    # Process new row
    print(f"Processing case {i+1} of {len(case_table)}")
    print(row)
    print()
    color = "black" if i == 0 else colors[i - 1]
    print(row["tag"])
    dir_out_current = os.path.join(DIR_OUT, row["tag"])

    # Create barotropic model object
    fluid_name = row["fluid_name"]
    model = bpy.BarotropicModel(
        fluid_name=fluid_name,
        T_in=row["T0_in"],
        p_in=row["p0_in"],
        p_out=row["p_out"],
        efficiency=row["polytropic_efficiency"],
        calculation_type=row["calculation_type"],
        # calculation_type="equilibrium",
        # blending_onset=row["q_onset"],
        # blending_width=row["q_transition"],
        blending_onset=0.0,
        blending_width=0.02,
        HEOS_solver="hybr",
        ODE_solver="LSODA",
        ODE_tolerance=1e-9,
        polynomial_degree=8,
        polynomial_format="horner",
        output_dir=dir_out_current,
    )

    # Evaluate barotropic model and export polynomial expressions
    model.solve()
    model.fit_polynomials()


#     model.export_expressions_fluent(output_dir=dir_out)
#     model.export_expressions_cfx(output_dir=dir_out)
    # for var in model.poly_fitter.variables:
    #     model.poly_fitter.plot_polynomial_and_error(var=var, savefig=False, showfig=True)

    model.poly_fitter.plot_polynomial_and_error(var="density", savefig=False, showfig=True)


    p = np.linspace(model.poly_fitter.p_out/2, model.poly_fitter.p_in, 10000)
    rho = model.poly_fitter.evaluate_polynomial(p, "density")
    slope = (p[1:] - p[0:-1])/(rho[1:] - rho[0:-1])
    speed_sound = np.sqrt(slope)


    # Plot phase diagram
    fluid = bpy.Fluid(name=fluid_name, backend="HEOS", exceptions=True)
    fig, ax = plt.subplots(figsize=(6, 5.0))
    ax.set_xlabel("Pressure")
    ax.set_ylabel("Density")


    ax.plot(
        p,
        rho,
        # linewidth=1.0,
        # marker="o",
        # markersize=3,
        color=color,
        linewidth=1.25,
        label=rf"$q_\text{{onset}}={row['q_onset']:0.2f}$, {row['calculation_type']}",
    )

    # Plot phase diagram
    fluid = bpy.Fluid(name=fluid_name, backend="HEOS", exceptions=True)
    fig, ax = plt.subplots(figsize=(6, 5.0))
    ax.set_xlabel("Pressure")
    ax.set_ylabel("Speed of sound")

    ax.plot(
        p[0:-1],
        speed_sound,
        # linewidth=1.0,
        # marker="o",
        # markersize=3,
        color=color,
        linewidth=1.25,
        label=rf"$q_\text{{onset}}={row['q_onset']:0.2f}$, {row['calculation_type']}",
    )


    model.poly_fitter.plot_phase_diagram(
        fluid=fluid,
        var_x="s",
        var_y="T",
        output_dir=dir_out_current)

# Show figure
plt.show()
