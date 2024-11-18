import os
import shutil
import datetime
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import ansys.fluent.core as pyfluent

import barotropy as bpy

bpy.print_package_info()
bpy.set_plot_options()

# Main directory
MAIN_DIR = "output"
if not os.path.exists(MAIN_DIR):
    os.makedirs(MAIN_DIR)

# Barotropic figures settings
SHOW_FIG = False
SAVE_FIG = True

# Fluent Parameters
N_ITER = 500  # Number of iterations
PROCESSORS_NUMBER = 16  # Number of processors for fluent
# PROCESSORS_NUMBER = os.cpu_count() - 2  # Use all cores except 2
FLUENT_FILE = "elliot_1968.cas.h5"
EXCEL_DATAFILE = "elliot_1968_data.xlsx"  # Case summary file
PLOT_REALTIME_RESIDUALS = True  # plot the residuals in real time in a Python figure

# Integration timescale
TIME_SCALE_FACTOR = 0.75

# Residual
tol = 1e-6
RES_CONTINUITY = tol
RES_X_VELOCITY = tol
RES_Y_VELOCITY = tol
RES_K = tol
RES_OMEGA = tol

# Relaxation factors
RELF_DENSITY = 0.05
RELF_K = 0.75
RELF_OMEGA = 0.75
RELF_BODYFORCE = 1.0
RELF_TURB_VISCOSITY = 1.0

# Turbulence boundary conditions
INLET_TURBULENT_INTENSITY = 0.05
INLET_TURBULENT_VISCOSITY_RATIO = 10

# Load cases from Excel
data = pd.read_excel(EXCEL_DATAFILE, sheet_name="SI Units")
case_indices = [1, 2, 14]  # Consider only some of the cases
case_data = data[case_data["index"].isin(case_indices)]

# Opening Fluent
solver = pyfluent.launch_fluent(
    product_version=pyfluent.FluentVersion.v242,
    mode=pyfluent.FluentMode.SOLVER,
    ui_mode=pyfluent.UIMode.HIDDEN_GUI,
    dimension=pyfluent.Dimension.TWO,
    precision=pyfluent.Precision.DOUBLE,
    processor_count=PROCESSORS_NUMBER,
)


# Loop over cases
for i, row in case_data.iterrows():

    # Unique datetime for current case
    current_datetime = datetime.datetime.now()
    formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H-%M-%S")

    # Create case folder
    case_index = row["index"]
    case_name = f"case_{case_index}"
    case_dir = os.path.join(MAIN_DIR, case_name)
    os.makedirs(case_dir, exist_ok=True)

    # Timestamped subdirectory for the current simulation
    timestamp_dir = os.path.join(case_dir, formatted_datetime)
    os.makedirs(timestamp_dir, exist_ok=True)

    # ------------------------------ #
    # ------ Barotropic Model ------ #
    # ------------------------------ #
    print("\n")
    print(f"    Creating barotropic model for case {case_index}")
    print("\n")

    # Create Barotropic Model Subfolder
    barotropic_dir = os.path.join(timestamp_dir, f"barotropic_model")
    os.makedirs(barotropic_dir, exist_ok=True)

    # Create barotropic model object
    model = bpy.BarotropicModel(
        fluid_name=[row["fluid_1"], row["fluid_2"]],
        mixture_ratio=row["mixture_ratio"],
        T_in=row["T_in"],
        p_in=row["p_in"],
        p_out=row["p_out"],
        efficiency=row["polytropic_efficiency"],
        ODE_solver="LSODA",
        ODE_tolerance=1e-9,
        polynomial_degree=8,
        polynomial_format="horner",
        output_dir=barotropic_dir,
    )

    # Evaluate barotropic model and export polynomial expressions
    model.solve()
    model.fit_polynomials()
    model.export_expressions_fluent(output_dir=barotropic_dir)
    model.export_expressions_cfx(output_dir=barotropic_dir)
    for var in model.poly_fitter.variables:
        model.poly_fitter.plot_polynomial_and_error(
            var=var, savefig=SAVE_FIG, showfig=SHOW_FIG
        )

    # ------------------------------- #
    # ------ Fluent Simulation ------ #
    # ------------------------------- #

    # Create simulation_subfolder
    fluent_dir = os.path.join(timestamp_dir, f"fluent_simulation")
    os.makedirs(fluent_dir, exist_ok=True)

    print("\n")
    print(f"    Running Fluent simulations for case {case_index}")
    print("\n")

    transcript_file = os.path.join(fluent_dir, f"case_{case_index}.trn")
    solver.file.start_transcript(file_name=transcript_file)

    if PLOT_REALTIME_RESIDUALS:
        subprocess = bpy.plot_residuals_real_time(transcript_file)

    # Upload the case file
    solver.file.read_case(file_name=FLUENT_FILE)

    #  Inlet condition
    inlet = solver.setup.boundary_conditions.pressure_inlet["inlet"]
    inlet.momentum.supersonic_or_initial_gauge_pressure.value = row["p_in"]
    inlet.momentum.gauge_total_pressure.value = row["p_in"]
    inlet.turbulence.turbulent_intensity = INLET_TURBULENT_INTENSITY
    inlet.turbulence.turbulent_viscosity_ratio = INLET_TURBULENT_VISCOSITY_RATIO

    # Outlet condition
    outlet = solver.setup.boundary_conditions.pressure_outlet["outlet"]
    outlet.momentum.gauge_pressure.value = row["p_out"]

    # Viscious model
    viscous = solver.setup.models.viscous
    viscous.model = "k-omega"
    viscous.k_omega_model = "sst"

    # Read fluent expressions
    filename = os.path.join(barotropic_dir, "fluent_expressions.txt")
    expressions = bpy.read_expression_file(filename)

    # Upload the expressions
    expr = solver.setup.named_expressions
    expr["barotropic_density"] = {"definition": expressions["barotropic_density"]}
    expr["barotropic_viscosity"] = {"definition": expressions["barotropic_viscosity"]}

    # Set residuals
    residuals = solver.solution.monitor.residual.equations
    residuals["continuity"].absolute_criteria = RES_CONTINUITY
    residuals["x-velocity"].absolute_criteria = RES_X_VELOCITY
    residuals["y-velocity"].absolute_criteria = RES_Y_VELOCITY
    residuals["k"].absolute_criteria = RES_K
    residuals["omega"].absolute_criteria = RES_OMEGA
    # solver.tui.solve.set.advanced.retain_cell_residuals('yes')

    # Set time step method
    solver.tui.solve.set.pseudo_time_method.advanced_options()
    time_step_method = (
        solver.solution.run_calculation.pseudo_time_settings.time_step_method
    )
    time_step_method.time_step_method = "automatic"
    time_step_method.length_scale_methods = "conservative"
    time_step_method.time_step_size_scale_factor = TIME_SCALE_FACTOR

    # Set relaxation factors
    relax_factors = "/solve/set/pseudo-time-method/relaxation-factors"
    solver.execute_tui(rf"""{relax_factors}/density {RELF_DENSITY}""")
    solver.execute_tui(rf"""{relax_factors}/k {RELF_K}""")
    solver.execute_tui(rf"""{relax_factors}/omega {RELF_OMEGA}""")
    solver.execute_tui(rf"""{relax_factors}/body-force {RELF_BODYFORCE}""")
    solver.execute_tui(rf"""{relax_factors}/turb-viscosity {RELF_TURB_VISCOSITY}""")

    # # Define the out file
    # TODO: @Davoud: We need to fix this part, it was giving problems for the Elliot1968 case
    # out_file = os.path.join(dir_simulation, f"report-file-0_{formatted_datetime}.out")
    # solver.solution.monitor.report_files["report-file-0"] = {"file_name": out_file}
    # solver.solution.monitor.report_files["report-file-0"] = {
    #     "report_defs": ["mass_in", "mass-balance"]
    # }

    # TODO: @Davoud: write code to compute integral quantities for Elliot1968 case

    # Initialization
    solver.solution.initialization.initialization_type = "standard"
    solver.solution.initialization.standard_initialize()
    solver.solution.initialization.initialization_type = "hybrid"
    solver.solution.initialization.hybrid_initialize()

    # Run Simulation
    solver.solution.run_calculation.iterate(iter_count=N_ITER)

    # Create pressure profile solution file
    vars = ["pressure", "density", "velocity-magnitude"]
    for var in vars:

        output_var_path = f"case_{case_index}_{var}.xy"
        output_var_path = os.path.join(fluent_dir, output_var_path)
        output_var_path = '"' + output_var_path + '"'

        # Create pressure results file
        solver.tui.plot.plot(
            "yes",
            output_var_path,
            "yes",
            "no",
            "no",
            var,
            "no",
            "no",
            "x-coordinate",
            "wall_conv",
            "wall_throat",
            "wall_dive",
            "centerline_conv",
            "centerline_throat",
            "centerline_dive",
            "()",
        )

    # Create plots from CFD data (comparison between CFD and EXP to be done in a different script)
    bpy.plot_nozzle_xy_data(case_name=case_name, dir_xy=fluent_dir, data_vars=vars)

    # Save case data
    output_file_case = rf"case_{case_index}.cas.h5"
    case_file_path = os.path.join(fluent_dir, output_file_case)
    solver.file.write(file_name=case_file_path, file_type="case-data")
    solver.file.stop_transcript()

    # Plot residual history
    fullpath = os.path.join(fluent_dir, f"{case_name}_residuals")
    bpy.plot_residuals(transcript_file, savefig=True, fullpath=fullpath)

    # Update the 'latest' folder for this case with the most recent results
    latest_dir = os.path.join(case_dir, "latest")
    if os.path.isdir(latest_dir):
        shutil.rmtree(latest_dir)  # Remove previous latest folder for this case
    shutil.copytree(timestamp_dir, latest_dir)

    print("\n")
    print(f"    Simulation Done")
    print("\n")

    if PLOT_REALTIME_RESIDUALS:
        subprocess.kill()


# Exiting Fluent
solver.exit()

print("All cases simulated ")
