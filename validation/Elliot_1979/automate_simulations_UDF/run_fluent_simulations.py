import os
import shutil
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from barotropy import (
    set_plot_options,
    create_fluent_journal,
    run_fluent_simulation,
    plot_residuals,
    plot_nozzle_data,
)

# Import data about simulation cases
case_data = pd.read_excel("simulation_cases.xlsx", skiprows=0)
case_data = case_data[case_data["index"].isin([6])]

# Define static filenames and paths
SIMULATIONS_ROOTDIR = "output"
TEMPLATE_JOURNAL_FILENAME = "template_nozzle_elliot_1979.jou"
TEMPLATE_CASE_FILENAME = "template_nozzle_elliot_1979_UDF.cas.h5"

# Define several iteration stages to solve the problem
# Useful to specify a progressive tightening of relaxation factors
RESIDUAL_TOLERANCE = 1e-9
SOLUTION_STRATEGY = {
    "time_scale_factors": [0.25, 0.24, 0.23, 0.22, 0.20],
    "density_relaxation_factors": [1.0, 1e-3, 1e-6, 1e-9, 1e-12],
    "number_of_iterations": [500, 250, 250, 500, 500],
}

# SOLUTION_STRATEGY = {
#     "time_scale_factors": [0.50],
#     "density_relaxation_factors": [1.0],
#     "number_of_iterations": [100],
# }

# Define mapping betweet export filenames and Fluent variable names
EXPORT_VARS = {
    "pressure": "pressure",
    "density": "density",
    "viscosity": "viscosity-lam",
    "speed_sound": "sound-speed",
    "velocity": "velocity-magnitude",
    "y_plus": "y-plus",
    "barotropic_density": "expr:barotropic_density",
    "barotropic_speed_sound": "expr:barotropic_speed_sound",
    "barotropic_viscosity": "expr:barotropic_viscosity",
    # "barotropic_void_fraction": "expr:barotropic_void_fraction",
    # "barotropic_mass_fraction": "expr:barotropic_mass_fraction",
}

# Define which barotropic model variables should be used
BAROTROPIC_VARS = [var for var in EXPORT_VARS if var.startswith("barotropic_")]

# Define data from which zones should be exported to XY files
EXPORT_ZONES = ["wall", "axis"]

# Define run parameters
N_PROC = 10  # Number of parallel processes for the simulation
TIMEOUT = 30 * 60  # Maximum time per simulation in seconds

# Loop over selected cases
for _, row in case_data.iterrows():

    # Get the current case name from the Identifier and Case columns of Excel file
    index = int(row["index"])
    tag = row["tag"]
    case_name = f"case_{index:03d}" if pd.isna(tag) else f"case_{index:03d}_{tag}"
    print(f"Processing {case_name}")

    # Determine the path for the current case based on the root directory and case name
    case_path = os.path.join(SIMULATIONS_ROOTDIR, case_name)

    # Create a new results directory named with unique date-time suffix
    timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    results_dir_timestamp = os.path.join(case_path, f"simulation_{timestamp}")
    os.makedirs(results_dir_timestamp)
    # results_dir_timestamp = os.path.join(case_path, "simulation_latest")

    # Define the mapping between templated variables and actual values
    file_prefix = os.path.join(results_dir_timestamp, case_name)
    template_map = {
        "case_file": TEMPLATE_CASE_FILENAME,
        "transcript_file": file_prefix + ".trn",
        "report_file": file_prefix + ".out",
        "case_out_file": file_prefix + ".cas.h5",
        # "p0_in": row["p0_in"],
        # "p_out": row["p_out"],
        # "turbulence_intensity_in": row["turbulence_intensity_in"],
        # "viscosity_ratio_in": row["viscosity_ratio_in"],
        # "wall_roughness_constant": row["wall_roughness_constant"],
        # "wall_roughness_height": row["wall_roughness_height"],
        "residual_tol": RESIDUAL_TOLERANCE,
    }

    # Create journal file for the current case
    journal = create_fluent_journal(
        case_name,
        results_dir_timestamp,
        TEMPLATE_JOURNAL_FILENAME,
        template_map,
        barotropic_vars=BAROTROPIC_VARS,
        export_vars=EXPORT_VARS,
        export_zones=EXPORT_ZONES,
        solution_strategy=SOLUTION_STRATEGY,
    )

    # Write the completed journal to a file for use with Fluent
    journal_file = os.path.join(results_dir_timestamp, case_name + ".jou")
    with open(journal_file, "w") as file:
        file.write(journal)

    # Run Fluent Simulations
    run_fluent_simulation(
        journal,
        n_proc=N_PROC,
        timeout=TIMEOUT,
        print_transcript=True,
        plot_residuals=True,
        transcript_file=template_map["transcript_file"],
    )

    # Create the folder to save figures
    set_plot_options(grid=False)
    fig_dir = os.path.join(results_dir_timestamp, "figures")
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    # Plot residuals
    filename_res = os.path.join(fig_dir, f"{case_name}_residuals")
    plot_residuals(template_map["transcript_file"], savefig=True, fullpath=filename_res)

    # Plot nozzle data
    plot_nozzle_data(case_name, fig_dir, EXPORT_VARS.keys())

    # Make a copy of the latest results in a folder with a simpler name
    results_backup_dir = os.path.join(case_path, "simulation_latest")
    if os.path.exists(results_backup_dir):
        shutil.rmtree(results_backup_dir)

    # Copy the contents of the latest folder with a new folder name
    shutil.copytree(results_dir_timestamp, results_backup_dir)



