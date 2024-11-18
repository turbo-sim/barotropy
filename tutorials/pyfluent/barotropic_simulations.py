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

###################################################################################################################
# Select cases number
CASES = [72, 73, 69]

# Main directory
MAIN_DIR = "output"
if not os.path.exists(MAIN_DIR):
    os.makedirs(MAIN_DIR)

# Barotropic figures settings
SHOW_FIG = False
SAVE_FIG = True

# Fluent Parameters
N_ITER = 5000  # Number of iterations
SHOW_GRAPHICS = True  # Show fluent gui
PROCESSORS_NUMBER = 16  # Number of processors for fluent
FLUENT_FILE = "simoneau_hendricks_1979.cas.h5"
EXCEL_DATAFILE = "cases_summary.xlsx"  # Case summary file
PLOT_REALTIME_RESIDUALS = True  # plot the residuals in real time in a python figure

# Integration timescale
TIME_SCALE_FACTOR = 0.5

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

###################################################################################################################

# Load cases from Excel
data = pd.read_excel(EXCEL_DATAFILE)
case_data = data[data["Case"].isin(CASES)]

# Opening Fluent
solver = pyfluent.launch_fluent(
    product_version=pyfluent.FluentVersion.v242,
    mode=pyfluent.FluentMode.SOLVER,
    ui_mode=pyfluent.UIMode.HIDDEN_GUI,
    dimension=pyfluent.Dimension.TWO,
    precision=pyfluent.Precision.DOUBLE,
    processor_count=PROCESSORS_NUMBER,
)


# Retrieve the Fluent version
fluent_version = solver.get_fluent_version()
print(f"Current Fluent version: {fluent_version}")

# Loop over cases
for i, row in case_data.iterrows():

    # Unique datetime for current case
    current_datetime = datetime.datetime.now()
    formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H-%M-%S")

    # Case number simulated
    cas = row["Case"]

    # Create case folder
    dir_out = os.path.join(MAIN_DIR, f"case_{cas}")
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)

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
    print(f"    Creating barotropic model for case {cas}")
    print("\n")

    # Create Barotropic Model Subfolder
    dir_barotropic = os.path.join(dir_out, f"barotropic_model_{formatted_datetime}")
    if not os.path.isdir(dir_barotropic):
        os.makedirs(dir_barotropic)

    # Create fluid object
    fluid_name = row["fluidname"]
    fluid = bpy.Fluid(name=fluid_name, backend="HEOS")
    p_in = row["P_0_in"]  # total preessure in Pa
    p_out = row["P_out"]  # Static pressure in Pa
    T_in = row["T_0_in"]  # Total temperature in K
    polytropic_efficiency = row["polytropic_efficiency"]
    calculation_type = row["calculation_type"]
    q_onset = row["q_onset"]
    q_transition = row["q_transition"]

    # Create barotropic model object
    model = bpy.BarotropicModel(
        fluid_name=fluid_name,
        T_in=T_in,
        p_in=p_in,
        p_out=fluid.triple_point_liquid.p,  # Use triple pressure as minimum pressure
        efficiency=polytropic_efficiency,
        calculation_type=calculation_type,
        blending_onset=q_onset,
        blending_width=q_transition,
        HEOS_solver="hybr",
        ODE_solver="LSODA",
        ODE_tolerance=1e-9,
        polynomial_degree=8,
        polynomial_format="horner",
        output_dir=dir_barotropic,
        polynomial_variables=["density", "viscosity", "speed_sound"],
    )

    # Evaluate barotropic model and export polynomial expressions
    model.solve()
    model.fit_polynomials()
    model.export_expressions_fluent(output_dir=dir_barotropic)
    model.export_expressions_cfx(output_dir=dir_barotropic)
    model.poly_fitter.plot_phase_diagram(
        fluid=fluid, var_x="s", var_y="T", savefig=SAVE_FIG, showfig=SHOW_FIG
    )
    for var in model.poly_fitter.variables:
        model.poly_fitter.plot_polynomial_and_error(
            var=var, savefig=SAVE_FIG, showfig=SHOW_FIG
        )

    # ------------------------------- #
    # ------ Fluent Simulation ------ #
    # ------------------------------- #

    dir_simulation = os.path.join(dir_out, f"simulation_{formatted_datetime}")
    if not os.path.isdir(dir_simulation):
        os.makedirs(dir_simulation)

    print("\n")
    print(f"    Running Fluent simulations for case {cas}")
    print("\n")

    transcript_file = os.path.join(dir_simulation, f"case_{cas}.trn")
    solver.file.start_transcript(file_name=transcript_file)

    if PLOT_REALTIME_RESIDUALS:
        subprocess = bpy.plot_residuals_real_time(transcript_file)

    # Upload the case file
    # TODO: @Alberto: 01.11.2024 it seems that Fluent is changing the .cas file when the commands are being applied
    # TODO: I suggest to create a copy of the file in the folder corresponding to the current simulation case so that the commands are applied to the copy of the file
    solver.file.read_case(file_name=FLUENT_FILE)

    #  Inlet condition
    inlet = solver.setup.boundary_conditions.pressure_inlet["inlet"]
    inlet.momentum.supersonic_or_initial_gauge_pressure.value = p_in
    inlet.momentum.gauge_total_pressure.value = p_in
    inlet.turbulence.turbulent_intensity = row["inlet_turbulent_intensity"]
    inlet.turbulence.turbulent_viscosity_ratio = row["inlet_turbulent_viscosity_ratio"]

    # Outlet condition
    outlet = solver.setup.boundary_conditions.pressure_outlet["outlet"]
    outlet.momentum.gauge_pressure.value = p_out

    # Viscious model
    viscous = solver.setup.models.viscous
    viscous.model = "k-omega"
    viscous.k_omega_model = "sst"

    # Read fluent expressions
    filename = os.path.join(dir_barotropic, "fluent_expressions.txt")
    expressions = bpy.read_expression_file(filename)

    # Upload the expressions
    solver.setup.named_expressions["barotropic_density"] = {
        "definition": expressions["barotropic_density"]
    }
    solver.setup.named_expressions["barotropic_viscosity"] = {
        "definition": expressions["barotropic_viscosity"]
    }

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
    solver.solution.run_calculation.pseudo_time_settings.time_step_method.time_step_method = (
        "automatic"
    )
    solver.solution.run_calculation.pseudo_time_settings.time_step_method.length_scale_methods = (
        "conservative"
    )
    solver.solution.run_calculation.pseudo_time_settings.time_step_method.time_step_size_scale_factor = (
        TIME_SCALE_FACTOR
    )

    # Set relaxation factors
    solver.execute_tui(
        rf"""/solve/set/pseudo-time-method/relaxation-factors/density {RELF_DENSITY}"""
    )
    solver.execute_tui(
        rf"""/solve/set/pseudo-time-method/relaxation-factors/k {RELF_K}"""
    )
    solver.execute_tui(
        rf"""/solve/set/pseudo-time-method/relaxation-factors/omega {RELF_OMEGA}"""
    )
    solver.execute_tui(
        rf"""/solve/set/pseudo-time-method/relaxation-factors/body-force {RELF_BODYFORCE}"""
    )
    solver.execute_tui(
        rf"""/solve/set/pseudo-time-method/relaxation-factors/turb-viscosity {RELF_TURB_VISCOSITY}"""
    )

    # Define the out file
    out_file = os.path.join(dir_simulation, f"case_{cas}.out")
    solver.solution.monitor.report_files["report-file-0"] = {"file_name": out_file}
    solver.solution.monitor.report_files["report-file-0"] = {
        "report_defs": ["mass_in", "mass-balance"]
    }

    # Initialization
    solver.solution.initialization.initialization_type = "standard"
    solver.solution.initialization.standard_initialize()
    solver.solution.initialization.initialization_type = "hybrid"
    solver.solution.initialization.hybrid_initialize()

    # Run Simulation
    solver.solution.run_calculation.iterate(iter_count=N_ITER)

    # Create pressure profile solution file
    output_press_file = f"case_{cas}.xy"
    file_path_pressure = os.path.join(dir_simulation, output_press_file)
    output_press_file = '"' + file_path_pressure + '"'

    while os.path.exists(file_path_pressure):
        print(f"The file '{output_press_file}' exists.")

        # Pause the execution and ask for user input
        input("Delete the file and press Enter to continue.")

    else:
        # Create pressure results file
        solver.tui.plot.plot(
            "yes",
            output_press_file,
            "yes",
            "no",
            "no",
            "pressure",
            "no",
            "no",
            "x-coordinate",
            "wall",
            "()",
        )

    # Save case data
    output_file_case = rf"case_{cas}.cas.h5"
    case_file_path = os.path.join(dir_simulation, output_file_case)

    while os.path.exists(case_file_path):
        print(f"The file '{output_file_case}' exists.")

        # Pause the execution and ask for user input
        input("Delete the file and press Enter to continue.")

    else:
        # Save the case file
        solver.file.write(file_name=case_file_path, file_type="case")

    print("\n")
    print(f"    Simulation Done")
    print("\n")

    solver.file.stop_transcript()

    # Save residual plot
    filename = os.path.join(dir_simulation,f"case_{cas}_residuals")
    bpy.plot_residuals(transcript_file, savefig=True, fullpath=filename)

    if PLOT_REALTIME_RESIDUALS:
        subprocess.kill()


# Exiting Fluent
solver.exit()

print(" All cases simulated ")
