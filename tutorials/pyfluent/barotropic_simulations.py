import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ansys.fluent.core as pyfluent
import datetime
import time

import barotropy as bpy

bpy.print_package_info()
bpy.set_plot_options()

###################################################################################################################
# Select cases number
CASES = [68,70]

# Main directory
MAIN_DIR = "output"
if not os.path.exists(MAIN_DIR):
    os.makedirs(MAIN_DIR)

# Barotropic figures settings
SHOW_FIG = False
SAVE_FIG = True

# Fluent Parameters
N_ITER = 30 # Number of iterations
SHOW_GRAPHICS = True # Show fluent gui
PROCESSORS_NUMBER = 2 # Number of processors for fluent
FLUENT_FILE =  "simoneau_hendricks_1979.cas.h5"
EXCEL_DATAFILE = "cases_summary.xlsx"  # Case summary file
PLOT_REALTIME_RESIDUALS = False # plot the residuals in real time in a python figure

# Residual
RES_CONTINUITY = 1e-7
RES_X_VELOCITY = 1e-7
RES_Y_VELOCITY = 1e-7
RES_K = 1e-7
RES_OMEGA = 1e-7
TIME_SCALE_FACTOR = 0.05

# Relaxation factors
RELF_DENSITY = 0.01
RELF_K = 0.75
RELF_OMEGA = 0.75
RELF_BODYFORCE = 1.
RELF_TURB_VISCOSITY = 1.


###################################################################################################################

# Read Case Summary
data = pd.read_excel(EXCEL_DATAFILE)
case_data = data[data["Case"].isin(CASES)]

# Opening Fluent
solver = pyfluent.launch_fluent(show_gui = SHOW_GRAPHICS,   # to launch also the gui
                                precision = 'double',
                                processor_count = PROCESSORS_NUMBER,
                                version= '2d'
                                )





for i, row in case_data.iterrows():

    # Case number simulated
    cas = row["Case"]

    print("\n")
    print(f"    Case {cas}      ")
    print("\n")

    # Create case folder
    dir_out = os.path.join(MAIN_DIR, f"case_{cas}")
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    

    #####   Barotropic Model    #####

    # Create Barotropic Model Subfolder
    dir_barotropic = os.path.join(dir_out, f"barotropic_model")
    if not os.path.isdir(dir_barotropic):
        os.makedirs(dir_barotropic)
    

    # Create fluid object
    fluid_name = row["fluidname"]
    fluid = bpy.Fluid(name=fluid_name, backend="HEOS")

    # Define case parameters
    dT_subcooling = 10

    # total preessure in Pa
    p_in = row["P_0_in"]
    # Static pressure in Pa
    p_out = row["P_out"]
    # Total temperature in K
    T_in = row["T_0_in"]

    polytropic_efficiency = row["polytropic_efficiency"]
    calculation_type = row["calculation_type"]
    q_onset = row["q_onset"]
    q_transition = row["q_transition"]


    # Create barotropic model object
    model = bpy.BarotropicModel(
        fluid_name=fluid_name,
        T_in=T_in,
        p_in=p_in,
        p_out=p_out/2,
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
        polynomial_variables=['density', 'viscosity', 'speed_sound']
    )

    # Evaluate barotropic model and export polynomial expressions
    model.solve()
    model.fit_polynomials()
    model.export_expressions_fluent(output_dir=dir_barotropic)
    model.export_expressions_cfx(output_dir=dir_barotropic)
    for var in model.poly_fitter.variables:
        model.poly_fitter.plot_polynomial_and_error(var=var, savefig=SAVE_FIG, showfig=SHOW_FIG)

    # Plot phase diagram
    var_x = "s"
    var_y = "T"
    fig, (ax_1, ax_2) = plt.subplots(1, 2, figsize=(12.0, 5.0), gridspec_kw={"wspace": 0.25})
    fig.suptitle(f"Barotropic model for expansion of {fluid_name}", fontsize=14)
    ax_1.set_xlabel("Entropy (J/kg/K)\n")
    ax_1.set_ylabel("Temperature (K)")
    ax_2.set_xlabel("Pressure (Pa)\n")
    ax_2.set_ylabel("Density (kg/m$^3$)")
    ax_1 = fluid.plot_phase_diagram(
        var_x,
        var_y,
        axes=ax_1,
        plot_critical_point=True,
        plot_saturation_line=True,
        plot_spinodal_line=True,
        plot_quality_isolines=True,
        N=50,
    )

    ax_1.plot(
        model.states[var_x],
        model.states[var_y],
        color=bpy.COLORS_MATLAB[0],
        linewidth=1.25,
    )

    ax_2.plot(
        model.states["p"],
        model.states["density"],
        linewidth=1.25,
        marker="none",
        markersize=3,
        linestyle="-",
        color=bpy.COLORS_MATLAB[0],
        label="Calculated data"
    )

    p = np.linspace(p_out, p_in, 1000)
    
    ax_2.plot(
        p,
        model.poly_fitter.evaluate_polynomial(p, "density"),
        linewidth=1.25,
        marker="none",
        markersize=3,
        linestyle="-.",
        color=bpy.COLORS_MATLAB[1],
        label="Polynomial fit"
    )

    ax_2.legend(loc="lower right")

    # Save figure
    plt.savefig(os.path.join(dir_barotropic,"barotropic_model_expansion.png"))





    #####   Fluent Simulation   #####

    # Create Simulation Subfolder
    current_datetime = datetime.datetime.now()
    formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H-%M-%S")

    dir_simulation = os.path.join(dir_out, f"simulation_{formatted_datetime}")
    if not os.path.isdir(dir_simulation):
        os.makedirs(dir_simulation)
   

    print("\n")
    print(f"    Simulating Case {cas}")
    print("\n")

    transcript_file = os.path.join(dir_simulation,f"transcript_{formatted_datetime}.trn")
    solver.file.start_transcript(file_name = transcript_file)
    
    if PLOT_REALTIME_RESIDUALS:
        bpy.plot_residuals_real_time(transcript_file)

    # Upload the case file
    # TODO: @Alberto: 01.11.2024 it seems that Fluent is changing the .cas file when the commands are being applied
    # TODO: I suggest to create a copy of the file in the folder corresponding to the current simulation case so that the commands are applied to the copy of the file
    solver.file.read_case(file_name = FLUENT_FILE)

    #  Inlet condition
    inlet = solver.setup.boundary_conditions.pressure_inlet['inlet']
    inlet.momentum.supersonic_or_initial_gauge_pressure.value = p_in
    inlet.momentum.gauge_total_pressure.value = p_in
    inlet.turbulence.turbulent_intensity = row["inlet_turbulent_intensity"] 
    inlet.turbulence.turbulent_viscosity_ratio = row["inlet_turbulent_viscosity_ratio"]

    # Outlet condition
    outlet = solver.setup.boundary_conditions.pressure_outlet['outlet']
    outlet.momentum.gauge_pressure.value = p_out

    # Viscious model
    viscous = solver.setup.models.viscous
    viscous.model = "k-omega"
    viscous.k_omega_model = "sst"

    # Read fluent expressions
    filename = os.path.join(dir_barotropic,"fluent_expressions.txt")
    expressions = bpy.read_expression_file(filename)

    # Upload the expressions
    solver.setup.named_expressions['barotropic_density'] = {"definition" : expressions["barotropic_density"]}
    solver.setup.named_expressions['barotropic_viscosity'] = {"definition" : expressions["barotropic_viscosity"]}

    # Set residuals
    solver.solution.monitor.residual.equations['continuity'].absolute_criteria = RES_CONTINUITY
    solver.solution.monitor.residual.equations['x-velocity'].absolute_criteria = RES_X_VELOCITY
    solver.solution.monitor.residual.equations['y-velocity'].absolute_criteria = RES_Y_VELOCITY
    solver.solution.monitor.residual.equations['k'].absolute_criteria = RES_K
    solver.solution.monitor.residual.equations['omega'].absolute_criteria = RES_OMEGA
    # solver.tui.solve.set.advanced.retain_cell_residuals('yes')

    # Set time step method
    solver.tui.solve.set.pseudo_time_method.advanced_options()
    solver.solution.run_calculation.pseudo_time_settings.time_step_method.time_step_method = "automatic"
    solver.solution.run_calculation.pseudo_time_settings.time_step_method.length_scale_methods = "conservative"
    solver.solution.run_calculation.pseudo_time_settings.time_step_method.time_step_size_scale_factor = TIME_SCALE_FACTOR
    
    # Set relaxation factors
    solver.execute_tui(rf'''/solve/set/pseudo-time-method/relaxation-factors/density {RELF_DENSITY}''')
    solver.execute_tui(rf'''/solve/set/pseudo-time-method/relaxation-factors/k {RELF_K}''')
    solver.execute_tui(rf'''/solve/set/pseudo-time-method/relaxation-factors/omega {RELF_OMEGA}''')
    solver.execute_tui(rf'''/solve/set/pseudo-time-method/relaxation-factors/body-force {RELF_BODYFORCE}''')
    solver.execute_tui(rf'''/solve/set/pseudo-time-method/relaxation-factors/turb-viscosity {RELF_TURB_VISCOSITY}''')

    # Define the out file
    out_file = os.path.join(dir_simulation, f"report-file-0_{formatted_datetime}.out")
    solver.solution.monitor.report_files['report-file-0'] = {"file_name" : out_file}
    solver.solution.monitor.report_files['report-file-0'] = {"report_defs" : ["mass_in", "mass-balance"]}  

    # Initialization
    solver.solution.initialization.initialization_type = "standard"
    solver.solution.initialization.standard_initialize()
    solver.solution.initialization.initialization_type = "hybrid"
    solver.solution.initialization.hybrid_initialize()


    # Run Simulation
    solver.solution.run_calculation.iterate(iter_count = N_ITER)


    # Create pressure profile solution file
    output_press_file =  f'case_{cas}_{formatted_datetime}_pressure_profile.xy'
    file_path_pressure = os.path.join(dir_simulation,output_press_file)
    output_press_file = '"' + file_path_pressure + '"'


    while os.path.exists(file_path_pressure):
        print(f"The file '{output_press_file}' exists.")

        # Pause the execution and ask for user input
        input("Delete the file and press Enter to continue.")
    
    else:
        # Create pressure results file
        solver.tui.plot.plot('yes',
                        output_press_file,
                        'yes', 
                        'no', 
                        'no', 
                        'pressure', 
                        'no',
                        'no',
                        'x-coordinate',
                        'wall',
                        '()'
                         )

    # Save case data
    output_file_case =  rf'case_{cas}_{formatted_datetime}_simulation_setup.cas.h5'
    case_file_path = os.path.join(dir_simulation,output_file_case)

    while os.path.exists(case_file_path):
        print(f"The file '{output_file_case}' exists.")

        # Pause the execution and ask for user input
        input("Delete the file and press Enter to continue.")

    else:
         # Save the case file
        solver.file.write(file_name=case_file_path, 
                          file_type="case")


    print("\n")
    print(f"    Simulation Done")
    print("\n")

    solver.file.stop_transcript()

    if PLOT_REALTIME_RESIDUALS:
        plt.close('all')



# Exiting Fluent
solver.exit()

print(" All cases simulated ")