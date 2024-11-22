import os
import sys
import shutil
import signal
import datetime
import numpy as np
import pandas as pd
import pexpect.popen_spawn


from barotropy import (
    set_plot_options,
    render_fluent_journal,
    read_expression_file,
    plot_residuals_real_time,
)

# Define which cases to simulate
# case_indices = [90, 91, 92, 93, 94, 95]
case_indices = list(range(88, 96)) + list(range(100, 118))
case_indices = [100, 101]
case_indices = [92]

# Define static filenames and paths
JOURNALS_DIR = "./rendered_journals"
SIMULATIONS_DIR = "..\\data_simulations"
TEMPLATE_CASE_FILENAME = "template_simoneau_hendricks_1979.cas.h5"
TEMPLATE_JOURNAL_FILENAME = "template_simoneau_hendricks_1979.jou"

# Define turbulence quantities for the boundary conditions
TURBULENCE_INTENSITY = 5  # As percentage %
VISCOSITY_RATIO = 10  # Ratio of turbulence to molecular visocisities

# Define run parameters
N_PROC = 10  # Number of parallel processes for the simulation
MAX_ITER = 2000  # Maximum number of iterations
TIMEOUT = 20 * 60  # Maximum time per simulation in seconds

# Get path of Fluent executable
FLUENT_EXECUTABLE = os.getenv("FLUENT_BIN")

# Create a folder to save the rendered journal files if it does not exist
if not os.path.exists(JOURNALS_DIR):
    os.makedirs(JOURNALS_DIR)

# Get the case names (this directory should be populated by 'create_barotropic_models.m')
cases = os.listdir(SIMULATIONS_DIR)

# Loop over selected cases
for index in case_indices:
    # Find current case name using index as identifier
    case_index = f"{index:03d}"
    case_name = [case for case in cases if case.startswith(f"case_{index:03d}")].pop()
    case_path = os.path.join(SIMULATIONS_DIR, case_name)
    print(f"Processing {case_name}")

    # Create folder with a unique date-time identifies to save simulation results
    timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    results_dir = os.path.join(case_path, f"{timestamp}")
    os.makedirs(results_dir)
    # results_dir = os.path.join(case_path, f"simulation_results")

    # Get case boundary conditions from Excel file
    df = pd.read_excel("../cases_summary.xlsx", skiprows=0)
    row = df[df["filename"] == case_name + ".csv"]
    T_0_in = row["T_0_in"].iat[0]
    P_0_in = row["P_0_in"].iat[0]
    P_out = row["P_out"].iat[0]

    # # Define the mapping between templated variables and actual values
    # template_map = {
    #     "case_file": TEMPLATE_CASE_FILENAME,
    #     "transcript_file": os.path.join(results_dir, f"case_{case_index}.trn"),
    #     "report_file": os.path.join(results_dir, f"case_{case_index}.out"),
    #     "case_out_file": os.path.join(results_dir, f"case_{case_index}.cas.h5"),
    #     "pressure_file": os.path.join(results_dir, f"case_{case_index}_pressure.xy"),
    #     "turbulence_intensity_in": TURBULENCE_INTENSITY,
    #     "viscosity_ratio_in": VISCOSITY_RATIO,
    #     "P_0_in": P_0_in,
    #     "P_out": P_out,
    #     "iter": MAX_ITER,
    # }

    # Define the mapping between templated variables and actual values
    template_map = {
        "case_file": TEMPLATE_CASE_FILENAME,
        "transcript_file": os.path.join(results_dir, f"case_{case_index}.trn"),
        "report_file": os.path.join(results_dir, f"case_{case_index}.out"),
        "case_out_file": os.path.join(results_dir, f"case_{case_index}.cas.h5"),
        "pressure_file": os.path.join(results_dir, f"case_{case_index}_pressure.xy"),
        "density_file": os.path.join(results_dir, f"case_{case_index}_density.xy"),
        "velocity_file": os.path.join(results_dir, f"case_{case_index}_velocity.xy"),
        "barotropic_density_file": os.path.join(results_dir, f"case_{case_index}_barotropic_density.xy"),
        "barotropic_speed_sound_file": os.path.join(results_dir, f"case_{case_index}_barotropic_speed_sound.xy"),
        "barotropic_viscosity_file": os.path.join(results_dir, f"case_{case_index}_barotropic_viscosity.xy"),
        "barotropic_vol_frac_2_file": os.path.join(results_dir, f"case_{case_index}_barotropic_vol_frac_2.xy"),
        "turbulence_intensity_in": TURBULENCE_INTENSITY,
        "viscosity_ratio_in": VISCOSITY_RATIO,
        "P_0_in": P_0_in,
        "P_out": P_out,
        "iter": MAX_ITER,
    }

    # Read Barotropic Model expression file
    expression_file = os.path.join(
        case_path,
        "barotropic_model",
        f"case_{case_index}_fluent_expressions.txt",
    )
    barotropic_model_map = read_expression_file(expression_file)

    # Add the Barotropic Model exoressuibs to the template mapping dictionary
    template_map = {**template_map, **barotropic_model_map}

    # Render journal template and save a copy in a file
    journal = render_fluent_journal(TEMPLATE_JOURNAL_FILENAME, template_map)
    journal_file = os.path.join(JOURNALS_DIR, case_name + ".jou")
    with open(journal_file, "w") as file:
        file.write(journal)

    # Run simulations
    try:

        # Create Fluent child process
        # 2dpp -> 2D solver with double precision arithmetic
        # -g -> Do not start up the graphical user interface
        # -tN -> Run the simulation using N processor cores
        process = pexpect.popen_spawn.PopenSpawn(
            f'{os.getenv("FLUENT_EXEC")} 2ddp -g -t{N_PROC}'
        )

        # Send sequence of commands to Fluent CLI
        process.send(journal)

        # Plot residuals from transcript file in real time
        transcript_file = os.path.join(results_dir, f"case_{case_index}.trn")
        res_process = plot_residuals_real_time(transcript_file)  # Function waits until the transcript file exists

        # Wait for the child process to finish
        exit_flag = process.expect([pexpect.EOF, pexpect.TIMEOUT], timeout=TIMEOUT)
        # print(process.before.decode('utf-8'))

        # Check execution status
        if exit_flag == 0:
            # Fluent execution reached the End Of File (EOF)
            print("Simulation completed\n")

        elif exit_flag == 1:
            # Time out reached before Fluent execution was complete
            # Kill Fluent using the automatically-generated clean-up script
            print("Timeout before execution was complete. Closing fluent with clean-up script")
            for filename in os.listdir(os.getcwd()):
                if filename.startswith("cleanup-fluent"):
                    os.system(filename)

        # Make a copy of the latest results in a folder with a simpler name
        results_backup_dir = os.path.join(case_path, "simulation_results")
        if os.path.exists(results_backup_dir):
            shutil.rmtree(results_backup_dir)

        # Copy the contents of the latest folder with a new folder name
        shutil.copytree(results_dir, results_backup_dir)

    except KeyboardInterrupt:
        print("Execution terminated by user. Closing fluent with clean-up script")
        for filename in os.listdir(os.getcwd()):
            if filename.startswith("cleanup-fluent"):
                os.system(filename)
        break

    finally:
        # Ensure the plot residuals process is terminated when the main script finishes
        if res_process:
            res_process.terminate()
