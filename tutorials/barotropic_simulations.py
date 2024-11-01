import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ansys.fluent.core as pyfluent
import datetime

import barotropy as bpy

bpy.print_package_info()
bpy.set_plot_options()


###################################################################################################################
# Select cases number
cases = [68,70,69]

# Main directory
main_dir = "cfd_solutions"
if not os.path.exists(main_dir):
    os.makedirs(main_dir)

show_fig = False
save_fig = True

# Fluent Parameters
n_it = 10
show_graphics = True # Show fluent gui
processors_number = 22 

FLUENT_FILE =  "../projects/Simoneau_Hendricks_1979/automate_simulations/template_simoneau_hendricks_1979.cas.h5"


current_path = os.getcwd()
fullpath = os.path.abspath(FLUENT_FILE)
print(fullpath)


###################################################################################################################


# Read Case Summary
excel_datafile = "../projects/Simoneau_Hendricks_1979/cases_summary.xlsx"
data = pd.read_excel(excel_datafile)
case_data = data[data["Case"].isin(cases)]


# # for i in range(len(cases)):
for i, row in case_data.iterrows():

    cas = row["Case"]

    print("\n")
    print(f"    Case {cas}      ")
    print("\n")

    # Create folder to save results
    DIR_OUT = os.path.join(main_dir, f"case_{cas}")
    if not os.path.isdir(DIR_OUT):
        os.makedirs(DIR_OUT)

    ### Barotropic Model ###

    # Create fluid object
    fluid_name = row["fluidname"]
    fluid = bpy.Fluid(name=fluid_name, backend="HEOS")

    # Experimental results
    mass_exp = row["m_exp"] # Mass flow experiments kg/s
    press_exp_file = row["filename"] # file containing the pressure profile results

    # Define case parameters
    dT_subcooling = 10

    # total preessure in Pa
    p_in = row["P_0_in"]
    # Static pressure in Pa
    p_out = row["P_out"]
    # Total temperature in K
    T_in = row["T_0_in"]

    # TODO put variables in the excel file
    # TODO possibility to add new column with tag to the Excel file. Check example in Elliot_1982/automate_simulations/simulation_cases.xlsx
    polytropic_efficiency = 1.00
    calculation_type = "blending"
    q_onset = 0.10 # use the liquid line till quality of 5%
    q_transition = 0.05 # transient from metastable liquid to equilibrium zone


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
        output_dir=DIR_OUT,
        polynomial_variables=['density', 'viscosity', 'speed_sound']
    )

    # Evaluate barotropic model and export polynomial expressions
    model.solve()
    model.fit_polynomials()
    model.export_expressions_fluent(output_dir=DIR_OUT)
    model.export_expressions_cfx(output_dir=DIR_OUT)
    for var in model.poly_fitter.variables:
        model.poly_fitter.plot_polynomial_and_error(var=var, savefig=save_fig, showfig=show_fig)

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
    # TODO use os.path.join
    fig_path = DIR_OUT + "/" + f"barotropic_model_expansion.png"
    plt.savefig(fig_path)




    ###     Fluent Simulation      ###

    current_datetime = datetime.datetime.now()
    formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H-%M-%S")


    print("\n")
    print(f"    Simulating Case {cas}")
    print("\n")

    # Launch fluent
    solver = pyfluent.launch_fluent(show_gui = show_graphics,   # to launch also the gui
                                    #  mode = 'solver',
                                     precision = 'double',
                                     processor_count = processors_number,
                                     version= '2d'
                                     )




    # Upload the case file
    solver.file.read_case(file_name = FLUENT_FILE)


    #  Inlet condition
    inlet = solver.setup.boundary_conditions.pressure_inlet['inlet']
    inlet.momentum.supersonic_or_initial_gauge_pressure.value = p_in
    inlet.momentum.gauge_total_pressure.value = p_in
    inlet.turbulence.turbulent_intensity = 0.05  # TODO implement in excel file or at the beginning of the file
    inlet.turbulence.turbulent_viscosity_ratio = 10

    # Outlet condition
    outlet = solver.setup.boundary_conditions.pressure_outlet['outlet']
    outlet.momentum.gauge_pressure.value = p_out


    # Viscious model
    viscous = solver.setup.models.viscous
    viscous.model = "k-omega"
    viscous.k_omega_model = "sst"

    filename = DIR_OUT + "/fluent_expressions.txt" # TODO os.path.join

    # Read the barotropic functions
    # TODO implement as a loop to read expression fule only once and load all expressions
    # TODO, reimplement reader as function rather than class
    density = bpy.read_expression_file(filename=filename, var_name="density")
    density_expression = density.expression()

    viscosity = bpy.read_expression_file(filename=filename, var_name="viscosity")
    viscosity_expression = viscosity.expression()


    # Upload the expressions
    solver.setup.named_expressions['rhomass_barotropic'] = {"definition" : density_expression}
    solver.setup.named_expressions['viscosity_barotropic'] = {"definition" : viscosity_expression}


    # Initialization
    solver.solution.initialization.initialization_type = "standard"
    solver.solution.initialization.standard_initialize()
    solver.solution.initialization.initialization_type = "hybrid"
    solver.solution.initialization.hybrid_initialize()

    # Run Simulation
    solver.solution.run_calculation.iterate(iter_count = n_it)



    # Saving Folder

    save_fold = DIR_OUT + "/"

    # Create pressure profile solution file
    output_press2 =  rf'case_{cas}_{formatted_datetime}_pressure_profile.csv'
    output_press_file = '"' + save_fold + output_press2 + '"'
    file_path_pressure =  save_fold + output_press2

    # if os.path.exists(file_path_pressure):
    #     print(f"The file '{output_press2}' exists.")

    #     # Pause the execution and ask for user input
    #     input("Press Enter to continue or type 'exit' to quit: ")
    
    # else:
    #     # Create pressure results file
    #     solver.tui.plot.plot('yes',
    #                     output_press_file,
    #                     'yes', 
    #                     'no', 
    #                     'no', 
    #                     'pressure', 
    #                     'no',
    #                     'no',
    #                     'x-coordinate',
    #                     'wall',
    #                     '()'
    #                     )

    while os.path.exists(file_path_pressure):
        print(f"The file '{output_press2}' exists.")

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



    # Save massflow results
    # TODO: do not export flux calculations to .flp files
    # TODO: instead crete flux calculations with python code and add them to the .out file
    # TODO: when needed import the quantities of the .out file using the in-house function bpy.parse_fluent_out
    mass_output_file =  rf'case_{cas}_{formatted_datetime}_massflow_results.flp'
    file_path_mass =  save_fold + mass_output_file

    if os.path.exists(file_path_mass): ## TODO replace by while loop if needed
        print(f"The file '{mass_output_file}' exists.")
    
        # Pause the execution and ask for user input
        input("Press Enter to continue or type 'exit' to quit:")
    
    else:
        # Create mass flow results file
        solver.results.report.fluxes.mass_flow(zones = ["outlet"], 
                                           write_to_file = True,
                                           file_name = file_path_mass)
    

    # Compare Massflow
    flp_file = file_path_mass

    # Read the flp file
    try:
        with open(flp_file, 'r') as file:
            content = file.readlines()  # Read all lines into a list

        # Step 2: Extract the outlet mass flow rate
        outlet_mass_flow_rate = None  # Initialize variable to store the value
        for line in content:
            if "outlet" in line:  # Check if the line contains "outlet"
                # Split the line and extract the second element (the value)
                parts = line.split()
                if len(parts) >= 2:  # Ensure there are enough parts
                    outlet_mass_flow_rate = - float(parts[1])  # Convert to float
                    delta = outlet_mass_flow_rate - mass_exp
                    break  # Exit the loop after finding the value

        # Display the result
        if outlet_mass_flow_rate is not None:
            print(f"Outlet Mass Flow Rate: {outlet_mass_flow_rate} kg/s")
        else:
            print("Outlet Mass Flow Rate not found.")

    except FileNotFoundError:
        print(f"The file '{flp_file}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    # Rewrite the file adding the experimental data
    try:
        with open(flp_file, 'r') as file:
            content = file.readlines()  # Read all lines into a list

        # Check if there are at least 7 lines to replace
        if len(content) >= 7:
            # Step 2: Replace line 7 (index 6, since index starts from 0)
            content[6] = f"                      exp_mass       {mass_exp}\n"  # Replace with your desired text

            # Step 3: Add two new lines at the end
            content.append("-------------------------------- --------------------\n")
            content.append(f"          cfd_mass - exp_mass      {delta}\n")

            # Step 4: Write the modified content back to the file
            with open(flp_file, 'w') as file:
                file.writelines(content)  # Write all modified lines back to the file

        else:
            print("The file does not contain enough lines to replace line 7.")

    except FileNotFoundError:
        print(f"The file '{flp_file}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")



    # Save case data
    output_file =  rf'case_{cas}_{formatted_datetime}_simulation_setup.cas.h5'
    case_file = save_fold + output_file

    if os.path.exists(case_file):
        print(f"The file '{output_file}' exists.")
    
        # Pause the execution and ask for user input
        input("Press Enter to continue or type 'exit' to quit:")
    

    else:
         # Save the case file
        solver.file.write(file_name=case_file, 
                          file_type="case")



    solver.exit()

    print("\n")
    print(f"    Simulation Done")
    print("\n")

    ###### Compare Solution #######

    # Experimental Data
    exp_res = "../projects/Simoneau_Hendricks_1979/data_experiments/" + press_exp_file
    data = pd.read_csv(exp_res, skiprows=1, header=None)
    data.columns = ['Position', 'Pressure_Pa']
    plt.figure(figsize=(8, 6))
    plt.plot(data['Position'], data['Pressure_Pa'], marker='o', color='k', label='Experimental Data')

    #  Simulation data
    cfd_res = file_path_pressure
    data = bpy.parse_fluent_xy(file_path_pressure)  # TODO, make sure xy reader works wall, maybe save data at several surfaces like wall and centerline
    data = pd.read_csv(cfd_res, skiprows=4, sep='\s+', header=None)
    data.columns = ['Position', 'Pressure_Pa']
    data['Position'] = pd.to_numeric(data['Position'], errors='coerce')
    data['Pressure_Pa'] = pd.to_numeric(data['Pressure_Pa'], errors='coerce')
    plt.plot(data['Position'], data['Pressure_Pa'], linestyle='-', color='b', label='CFD Data')


    plt.xlabel('Position')
    plt.ylabel('Pressure (Pa)')
    plt.title(f"Case {cas} with barotropic model")
    plt.grid(True)
    plt.legend()


    fig_path = save_fold + f"case_{cas}_{formatted_datetime}_pressure_profile.png"
    plt.savefig(fig_path)
