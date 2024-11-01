
# %% Import Packages

import ansys.fluent.core as pyfluent
import pandas as pd
import matplotlib.pyplot as plt
import os
import datetime
import time
import barotropy as bpy

current_datetime = datetime.datetime.now()
formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H-%M-%S")

# Input Variables
cas = 73
n_it = 10
show_graphics = True


# Read Case Summary
case_summ = "../projects/Simoneau_Hendricks_1979/cases_summary.xlsx"
data = pd.read_excel(case_summ) # , skiprows=1)
data = pd.DataFrame(data)
case_data = data[data.iloc[:, 0] == cas]

# total preessure in Pa
p_in = float(case_data["P_0_in"])
# Static pressure in Pa
p_out = float(case_data["P_out"])
# Total temperature in K
T_in = float(case_data["T_0_in"])

# Experimental results
# Mass flow experiments kg/s
mass_exp = float(case_data["m_exp"])
# file containing the pressure profile results
press_exp_file = case_data["filename"].iloc[0] 

# Launch fluent
solver = pyfluent.launch_fluent(show_gui = show_graphics,   # to launch also the gui
                                #  mode = 'solver',
                                 precision = 'double',
                                 processor_count = 20,
                                 version= '2d'
                                 )




# Upload the case file
solver.file.read_case(file_name = "../projects/Simoneau_Hendricks_1979/automate_simulations/template_simoneau_hendricks_1979.cas.h5")



#  Inlet condition
inlet = solver.setup.boundary_conditions.pressure_inlet['inlet']
inlet.momentum.supersonic_or_initial_gauge_pressure.value = p_in
inlet.momentum.gauge_total_pressure.value = p_in
inlet.turbulence.turbulent_intensity = 0.05
inlet.turbulence.turbulent_viscosity_ratio = 10

# Outlet condition
outlet = solver.setup.boundary_conditions.pressure_outlet['outlet']
outlet.momentum.gauge_pressure.value = p_out


# Viscious model
viscous = solver.setup.models.viscous
viscous.model = "k-omega"
viscous.k_omega_model = "sst"

# Read the barotropic functions
# density = bpy.read_expression_file(filename="output/fluent_expressions.txt", var_name="density")
# density_expression = density.expression()

# viscosity = bpy.read_expression_file(filename="output/fluent_expressions.txt", var_name="viscosity")
# viscosity_expression = viscosity.expression()


expressions = bpy.read_expression_file(filename="output/fluent_expressions.txt")  # Function that returns a list of dictionaries
for expression in expressions:
    solver.setup.named_expressions[expression["description"]] = {"definition" : expression["value"]}


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

save_fold = "cfd_tests/"

# Create pressure profile solution file
output_press2 =  rf'case_{cas}_{formatted_datetime}_pressure_profile.csv'
output_press_file = '"' + save_fold + output_press2 + '"'
file_path_pressure =  save_fold + output_press2

if os.path.exists(file_path_pressure):
    print(f"The file '{output_press2}' exists.")
    
    # Pause the execution and ask for user input
    input("Press Enter to continue or type 'exit' to quit: ")
    
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
mass_output_file =  rf'case_{cas}_{formatted_datetime}_massflow_results.flp'
file_path_mass =  save_fold + mass_output_file

if os.path.exists(file_path_mass):
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

###### Compare Solution #######

# Experimental Data
exp_res = "../projects/Simoneau_Hendricks_1979/data_experiments/" + press_exp_file
data = pd.read_csv(exp_res, skiprows=1, header=None)
data.columns = ['Position', 'Pressure_Pa']
plt.figure(figsize=(8, 6))
plt.plot(data['Position'], data['Pressure_Pa'], marker='o', color='k', label='Experimental Data')

#  Simulation data
cfd_res = file_path_pressure
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

# plt.show()

# # Compare Massflow
# flp_file = path2 + mass_output_file  # Change to your .flp file name
# new_file = path2 + f"case_{cas}_{formatted_datetime}_massflow_comparison.txt"  # File to write modified content

# # Read the .flp file
# try:
#     with open(flp_file, 'r') as file:
#         content = file.readlines()  # Read all lines into a list

#     # Step 2: Extract the outlet mass flow rate
#     outlet_mass_flow_rate = None  # Initialize variable to store the value
#     for line in content:
#         if "outlet" in line:  # Check if the line contains "outlet"
#             # Split the line and extract the second element (the value)
#             parts = line.split()
#             if len(parts) >= 2:  # Ensure there are enough parts
#                 outlet_mass_flow_rate = - float(parts[1])  # Convert to float
#                 delta = outlet_mass_flow_rate - mass_exp
#                 break  # Exit the loop after finding the value

#     # Display the result
#     if outlet_mass_flow_rate is not None:
#         print(f"Outlet Mass Flow Rate: {outlet_mass_flow_rate} kg/s")
#     else:
#         print("Outlet Mass Flow Rate not found.")

# except FileNotFoundError:
#     print(f"The file '{flp_file}' was not found.")
# except Exception as e:
#     print(f"An error occurred: {e}")



# try:
#     with open(flp_file, 'r') as file:
#         content = file.readlines()  # Read all lines into a list

#     # Check if there are at least 7 lines to replace
#     if len(content) >= 7:
#         # Step 2: Replace line 7 (index 6, since index starts from 0)
#         content[6] = f"                      exp_mass       {mass_exp}\n"  # Replace with your desired text

#         # Step 3: Add two new lines at the end
#         content.append("-------------------------------- --------------------\n")
#         content.append(f"          cfd_mass - exp_mass      {delta}\n")

#         # Step 4: Write the modified content back to the file
#         with open(flp_file, 'w') as file:
#             file.writelines(content)  # Write all modified lines back to the file

#     else:
#         print("The file does not contain enough lines to replace line 7.")

# except FileNotFoundError:
#     print(f"The file '{flp_file}' was not found.")
# except Exception as e:
#     print(f"An error occurred: {e}")


