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
SHOW_FIG = True
SAVE_FIG = True

# Number of parallel processors for fluent
PROCESSORS_NUMBER = 24  
# PROCESSORS_NUMBER = os.cpu_count() - 2  # Use all cores except 2

# Residual tolerances
tol = 1e-4
RES_CONTINUITY = tol
RES_X_VELOCITY = tol
RES_Y_VELOCITY = tol
RES_Z_VELOCITY = tol
RES_K = tol
RES_OMEGA = tol
PLOT_REALTIME_RESIDUALS = True  # plot the residuals in real time in a Python figure

# Number of iterations
N_ITER = 2000  

# Integration timescale
TIME_SCALE_FACTOR = 5

# Relaxation factors
RELF_K = 0.75
RELF_OMEGA = 0.75
RELF_BODYFORCE = 1.0
RELF_TURB_VISCOSITY = 1.0
RELF_DENSITY = 0.2

# Define solver strategies with lists
# SOLUTION_STRATEGY = {
#     "N_ITER": [100, 200, 200, 250, 250, 250],
#     "TIME_SCALE_FACTOR": [1.0, 2.0, 3.0, 1.0, 0.5, 0.1],
#     "RELF_DENSITY": [0.10, 0.25, 0.50, 0.10, 0.05, 0.01],
# }
SOLUTION_STRATEGY = {
    "N_ITER": [200, 200, 200, 200, 200],
    "TIME_SCALE_FACTOR": [2.0, 1.5, 1.00, 0.75, 0.50],
    "RELF_DENSITY": [0.50, 0.25, 0.10, 0.05, 0.01],
}

# Create fluid
fluid_name = "CO2"
fluid = bpy.Fluid(name=fluid_name, backend="HEOS")

# Read case summary
DATAFILE = "./nakagawa_2009_cases.xlsx"
case_data = pd.read_excel(DATAFILE)

# Select cases of interest
# GEOMETRY = ["nozzle_2"]
# case_data = case_data[case_data["nozzle"].isin(GEOMETRY)]
# case_data = case_data[case_data["tag"] == "efficiency_sens"]
# case_data = case_data[case_data["tag"] == "roughness_sens"]
# case_data = case_data[(case_data["tag"] == "equilibrium_sens_2") | (case_data["tag"] == "equilibrium_sens_3")]
# case_data = case_data[np.isclose(case_data["p0_in"], 91e5)]
# cases = ["inviscid"]
# cases = [300, 301, 302, 303]
# cases = [400, 401, 402, 403]
# cases = [606, 607, 608, 609, 610]
# case_data = case_data[case_data["case"].isin(cases)]
# case_data = case_data[case_data["tag"] == "metastable_sens"]

cases = [2002, 2010, 2011, 2012]
case_data = case_data[case_data["case"].isin(cases)]

# Loop over cases
print("Cases scheduled for simulation:")
print(", ".join(case_data["case"].astype(str)))
solver = None
for i, row in case_data.iterrows():

    # Launch Fluent for the first case only
    if solver is None:
        print("Launching Fluent...")
        solver = pyfluent.launch_fluent(
            product_version=pyfluent.FluentVersion.v242,
            mode=pyfluent.FluentMode.SOLVER,
            ui_mode=pyfluent.UIMode.HIDDEN_GUI,
            dimension=pyfluent.Dimension.THREE,
            precision=pyfluent.Precision.DOUBLE,
            processor_count=PROCESSORS_NUMBER,
        )
        print("Fluent launched successfully.")

    # Unique datetime for current case
    current_datetime = datetime.datetime.now()
    formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H-%M-%S")

    # Create case folder
    index = row["case"]
    case_name = f"case_{index}"
    case_dir = os.path.join(MAIN_DIR, case_name)
    os.makedirs(case_dir, exist_ok=True)

    # Timestamped subdirectory for the current simulation
    timestamp_dir = os.path.join(case_dir, formatted_datetime)
    os.makedirs(timestamp_dir, exist_ok=True)

    # ------------------------------ #
    # ------ Barotropic Model ------ #
    # ------------------------------ #
    print("\n")
    print(f"    Creating barotropic model for case {index}")
    print("\n")

    # Create Barotropic Model Subfolder
    barotropic_dir = os.path.join(timestamp_dir, f"barotropic_model")
    os.makedirs(barotropic_dir, exist_ok=True)

    # Create barotropic model object
    model = bpy.BarotropicModel(
        fluid_name=fluid_name,
        T_in=row["T0_in"] + 273.15,
        p_in=row["p0_in"],
        p_out=1.05 * fluid.triple_point_liquid.p,
        efficiency=row["polytropic_efficiency"],
        calculation_type=row["calculation_type"],
        blending_onset=row["q_onset"],
        blending_width=row["q_transition"],
        ODE_solver="LSODA",
        ODE_tolerance=1e-12,
        polynomial_degree=[4, 6, 8],
        polynomial_format="horner",
        output_dir=barotropic_dir,
    )

    # Evaluate barotropic model and export polynomial expressions
    model.solve()
    model.fit_polynomials()
    model.export_expressions_fluent(output_dir=barotropic_dir)
    model.export_expressions_cfx(output_dir=barotropic_dir)
    for name in model.poly_fitter.variables:
        model.poly_fitter.plot_polynomial_and_error(
            var=name, savefig=SAVE_FIG, showfig=SHOW_FIG
        )

    # ------------------------------- #
    # ------ Fluent Simulation ------ #
    # ------------------------------- #

    # Create simulation_subfolder
    fluent_dir = os.path.join(timestamp_dir, f"fluent_simulation")
    os.makedirs(fluent_dir, exist_ok=True)

    print("\n")
    print(f"    Running Fluent simulations for case {index}")
    print("\n")

    transcript_file = os.path.join(fluent_dir, f"{case_name}.trn")
    solver.file.start_transcript(file_name=transcript_file)

    if PLOT_REALTIME_RESIDUALS:
        subprocess = bpy.plot_residuals_real_time(transcript_file)

    # Upload the case file
    solver.file.read_case(file_name=f"cases/{row["case_file"]}")

    # Define viscous model
    viscous_model = row["viscous_model"]
    solver.setup.models.viscous.model = viscous_model
    if viscous_model == "k-omega":
        solver.setup.models.viscous.k_omega_model = "sst"

    #  Inlet condition
    inlet_zones = [
        "velocity-inlet-28",
        "velocity-inlet-36",
        "velocity-inlet-42",
        "velocity-inlet-48",
        "velocity-inlet-52",
    ]
    for zone in inlet_zones:
        inlet = solver.setup.boundary_conditions.pressure_inlet[zone]
        inlet.momentum.supersonic_or_initial_gauge_pressure.value = row["p0_in"]
        inlet.momentum.gauge_total_pressure.value = row["p0_in"]
        if viscous_model != "inviscid":
            inlet.turbulence.turbulent_intensity = row["turbulence_intensity_in"]
            inlet.turbulence.turbulent_viscosity_ratio = row["viscosity_ratio_in"]

    # Outlet condition
    outlet_zones = [
        "outflow-121",
        "outflow-128",
        "outflow-133",
        "outflow-138",
        "outflow-141",
    ]
    for zone in outlet_zones:
        outlet = solver.setup.boundary_conditions.pressure_outlet[zone]
        outlet.momentum.gauge_pressure.value = row["p_out"]

    wall_zones = [
        "wall-35",
        "wall-41",
        "wall-47",
        "wall-51",
        "wall-61",
        "wall-66",
        "wall-71",
        "wall-74",
        "wall-83",
        "wall-88",
        "wall-93",
        "wall-96",
        "wall-105",
        "wall-110",
        "wall-115",
        "wall-118",
        "wall-127",
        "wall-132",
        "wall-137",
        "wall-140",
    ]
    for zone in wall_zones:
        if viscous_model != "inviscid":
            height_value = row["wall_roughness_height"]
            height_temp = {"turbulence": {"roughness_height": {"value": height_value}}}
            solver.setup.boundary_conditions.wall[zone] = height_temp

    # Read fluent expressions
    filename = os.path.join(barotropic_dir, "fluent_expressions.txt")
    expressions = bpy.read_expression_file(filename)

    # Upload the expressions
    expr = solver.setup.named_expressions
    for key, value in expressions.items():
        expr[key] = {"definition": value}

    # Set residuals
    residuals = solver.solution.monitor.residual.equations
    residuals["continuity"].absolute_criteria = RES_CONTINUITY
    residuals["x-velocity"].absolute_criteria = RES_X_VELOCITY
    residuals["y-velocity"].absolute_criteria = RES_Y_VELOCITY
    residuals["z-velocity"].absolute_criteria = RES_Z_VELOCITY
    if viscous_model != "inviscid":
        residuals["k"].absolute_criteria = RES_K
        residuals["omega"].absolute_criteria = RES_OMEGA

    # Keep residuals to inspect the flow in case of non-convergence
    solver.tui.solve.set.advanced.retain_cell_residuals("yes")

    # Set time step method
    solver.tui.solve.set.pseudo_time_method.advanced_options()
    time_step_method = (
        solver.solution.run_calculation.pseudo_time_settings.time_step_method
    )
    time_step_method.time_step_method = "automatic"
    time_step_method.length_scale_methods = "conservative"
    time_step_method.time_step_size_scale_factor = TIME_SCALE_FACTOR
    if viscous_model == "inviscid":
        time_step_method.time_step_size_scale_factor = 2 * TIME_SCALE_FACTOR    

    # Set relaxation factors
    relax_factors = "/solve/set/pseudo-time-method/relaxation-factors"
    solver.execute_tui(rf"""{relax_factors}/density {RELF_DENSITY}""")
    solver.execute_tui(rf"""{relax_factors}/body-force {RELF_BODYFORCE}""")
    if viscous_model == "k-omega":
        solver.execute_tui(rf"""{relax_factors}/k {RELF_K}""")
        solver.execute_tui(rf"""{relax_factors}/omega {RELF_OMEGA}""")
        solver.execute_tui(rf"""{relax_factors}/turb-viscosity {RELF_TURB_VISCOSITY}""")

    # Define the out file
    out_file = os.path.join(fluent_dir, f"{case_name}.out")
    solver.solution.monitor.report_files["report_file"] = {"file_name": out_file}

    # Initialization
    solver.solution.initialization.initialization_type = "standard"
    solver.solution.initialization.standard_initialize()
    solver.solution.initialization.initialization_type = "hybrid"
    solver.solution.initialization.hybrid_initialize()

    # Run simulation
    for step_idx in range(len(SOLUTION_STRATEGY["N_ITER"])):
        print(f"Running step {step_idx + 1}...")
        time_step_factor = SOLUTION_STRATEGY["TIME_SCALE_FACTOR"][step_idx]
        density_relax_fractor = SOLUTION_STRATEGY["RELF_DENSITY"][step_idx]
        max_iter = SOLUTION_STRATEGY["N_ITER"][step_idx]

        # Set parameters for the current step
        time_step_method.time_step_size_scale_factor = time_step_factor
        solver.execute_tui(rf"""{relax_factors}/density {density_relax_fractor}""")
        
        # Run simulation
        solver.solution.run_calculation.iterate(iter_count=max_iter)
        

    # Export xy files with flow solution
    vars = bpy.EXPORT_VARS
    if viscous_model == "inviscid":
        vars.pop("viscosity", None)
        vars.pop("y_plus", None)

    for name, fluent_name in vars.items():

        # File cannot be exported if it contains semicolon
        output_var_path = f"{case_name}_{name}.xy"
        output_var_path = os.path.join(fluent_dir, output_var_path)
        output_var_path = '"' + output_var_path + '"'

        # Create pressure results file
        solver.tui.plot.plot(
            "yes",
            output_var_path,
            "yes",
            "no",
            "no",
            fluent_name,
            "no",
            "no",
            "z-coordinate",
            "wall_centerline",
            "nozzle_centerline",
            "()",
        )

    # Create plots from CFD data (comparison between CFD and EXP to be done in a different script)
    bpy.plot_nozzle_xy_data(
        case_name=case_name,
        dir_xy=fluent_dir,
        data_vars=list(vars.keys()),
        zone_list=["wall_centerline"],
        x_var="Z-Coordinate",
    )

    # Save case data
    output_file_case = rf"{case_name}.cas.h5"
    case_file_path = os.path.join(fluent_dir, output_file_case)
    solver.file.write(file_name=case_file_path, file_type="case-data")
    solver.file.stop_transcript()

    # Plot residual history
    fullpath = os.path.join(fluent_dir, f"{case_name}_residuals")
    bpy.plot_residuals(transcript_file, savefig=True, fullpath=fullpath)

    # Update the 'latest' folder for this case with the most recent results
    latest_dir = os.path.join(case_dir, "latest")
    if os.path.isdir(latest_dir):  # Remove previous latest folder for this case
        shutil.rmtree(latest_dir, ignore_errors=True)
    shutil.copytree(timestamp_dir, latest_dir, dirs_exist_ok=True)

    print("\n")
    print(f"    Simulation Done")
    print("\n")

    if PLOT_REALTIME_RESIDUALS:
        subprocess.kill()

    print(f"Completed simulation for case {row['case']}.")

# Exiting Fluent
solver.exit()
print("All scheduled cases completed.")

