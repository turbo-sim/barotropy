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
    compare_contents_or_files
)

# Import data about simulation cases
case_data = pd.read_excel("simulation_cases.xlsx", skiprows=0)
case_data = case_data[case_data["index"].isin([1])]
case_data = case_data[case_data["tag"].isin(["a", "b", "c"])]

# Define static filenames and paths
OUTPUT_ROOTDIR = "output"
TEMPLATE_JOURNAL_FILENAME = "template_nozzle_elliot_1979.jou"
TEMPLATE_CASE_FILENAME = "template_nozzle_elliot_1979_withBL.cas.h5"

# Define several iteration stages to solve the problem
# Useful to specify a progressive tightening of relaxation factors
RESIDUAL_TOLERANCE = 1e-9
SOLUTION_STRATEGY = {
    "time_scale_factors": [0.25, 0.24, 0.23, 0.22, 0.20],
    "density_relaxation_factors": [1.0, 1e-3, 1e-6, 1e-9, 1e-12],
    "number_of_iterations": [500, 250, 250, 500, 1000],
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
    "barotropic_void_fraction": "expr:barotropic_void_fraction",
    "barotropic_mass_fraction": "expr:barotropic_mass_fraction",
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
    case_path = os.path.join(OUTPUT_ROOTDIR, case_name)

    # Check if the folder for the current case exists
    if not os.path.exists(case_path):
        # TODO: Logic to create the barotropic model is still pending
        raise FileNotFoundError(f"Folder {case_path} does not exist.")

    # Create a directory to save simulation results
    outdir_latest = os.path.join(case_path, "simulation_latest")
    if not os.path.exists(outdir_latest):
        os.makedirs(outdir_latest)

    # Define the mapping between templated variables and actual values
    file_prefix = os.path.join(outdir_latest, case_name)
    template_map = {
        "case_file": TEMPLATE_CASE_FILENAME,
        "transcript_file": file_prefix + ".trn",
        "report_file": file_prefix + ".out",
        "case_out_file": file_prefix + ".cas.h5",
        "p0_in": row["p0_in"],
        "p_out": row["p_out"],
        "turbulence_intensity_in": row["turbulence_intensity_in"],
        "viscosity_ratio_in": row["viscosity_ratio_in"],
        "wall_roughness_constant": row["wall_roughness_constant"],
        "wall_roughness_height": row["wall_roughness_height"],
        "residual_tol": RESIDUAL_TOLERANCE,
    }

    # Create journal file for the current case
    journal = create_fluent_journal(
        case_name,
        outdir_latest,
        TEMPLATE_JOURNAL_FILENAME,
        template_map,
        barotropic_vars=BAROTROPIC_VARS,
        export_vars=EXPORT_VARS,
        export_zones=EXPORT_ZONES,
        solution_strategy=SOLUTION_STRATEGY,
    )

    # Run Fluent simulations if the journal file has changed
    latest_journal_file = os.path.join(outdir_latest, case_name + ".jou")
    if compare_contents_or_files(journal, latest_journal_file):
        print(f"Detected identical settings for {case_name}. Skipping simulation...")
    else:

        # Write the completed journal to a file for use with Fluent
        shutil.rmtree(outdir_latest)  # Clean previous simulation
        os.makedirs(outdir_latest)
        with open(latest_journal_file, "w") as file:
            file.write(journal)

        # Define name of new directory to save copy of simulation results
        # Name based on unique timestamp identifies generated before the simulation
        timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        outdir_timestamp = os.path.join(case_path, f"simulation_{timestamp}")

        # Run fluent simulations using the current jounral file
        print(f"Running simulation for {case_name}")
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
        fig_dir = os.path.join(outdir_latest, "figures")
        if not os.path.exists(fig_dir):
            os.makedirs(fig_dir)

        # Plot residuals
        fig_name = os.path.join(fig_dir, f"{case_name}_residuals")
        plot_residuals(template_map["transcript_file"], savefig=True, fullpath=fig_name)

        # Plot nozzle data
        plot_nozzle_data(case_name, fig_dir, EXPORT_VARS.keys())

        # Copy the simulation results in a folder with a unique datetime identifier
        shutil.copytree(outdir_latest, outdir_timestamp)

        # Close figures to free memory
        plt.close('all')
