import os
import pandas as pd
import matplotlib.pyplot as plt

import barotropy as bpy

# Define plot settings
bpy.print_package_info()
bpy.set_plot_options()
SAVE_FIGURES = True
SHOW_FIGURES = True

# Create folder to save results
DIR_OUT = "./output"
if not os.path.isdir(DIR_OUT):
    os.makedirs(DIR_OUT)

# Load cases from Excel
case_table = pd.read_excel("./simulation_cases.xlsx")
case_indices = case_table["index"]
case_indices = [1]
case_table = case_table[case_table["index"].isin(case_indices)]

# Loop over cases
for i, (idx, row) in enumerate(case_table.iterrows()):

    # Process new row
    print(f"Processing case {i+1} of {len(case_table)}")
    print(row)
    print()
    dir_out = os.path.join(DIR_OUT, f"case_{row['index']}")

    # Create barotropic model object
    model = bpy.BarotropicModel(
        fluid_name=[row["fluid_1"], row["fluid_2"]],
        mixture_ratio=row["mixture_ratio"],
        T_in=row["T0_in"],
        p_in=row["p0_in"],
        p_out=row["p_out"],
        efficiency=row["efficiency"],
        ODE_solver="LSODA",
        ODE_tolerance=1e-9,
        polynomial_degree=8,
        polynomial_format="horner",
        output_dir=dir_out,
    )

    # Evaluate barotropic model and export polynomial expressions
    model.solve()
    model.fit_polynomials()
    model.export_expressions_fluent(output_dir=dir_out)
    model.export_expressions_cfx(output_dir=dir_out)
    for var in model.poly_fitter.variables:
        model.poly_fitter.plot_polynomial_and_error(var=var, savefig=SAVE_FIGURES, showfig=SHOW_FIGURES)

# Show figure
plt.show()
