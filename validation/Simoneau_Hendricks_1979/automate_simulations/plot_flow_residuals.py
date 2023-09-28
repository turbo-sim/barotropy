import os
from parse import parse
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd


from barotropy import (
    parse_fluent_xy,
    set_plot_options,
    parse_fluent_out,
    plot_residuals,
)


# Define plot settings
show_figures = True
set_plot_options(fontsize=14, grid=True, margin=0.00)

# Define directory paths
SIMULATIONS_DIR = "../data_simulations/"
EXPERIMENTAL_DIR = "../data_experiments/"
RESULTS_DIR = "simulation_results"

# Define which cases to plot
# case_indices = [90, 91, 92, 93, 94, 95]
case_indices = list(range(88, 96)) + list(range(100, 118))
case_indices = [112]

# Get the case names
cases = os.listdir(SIMULATIONS_DIR)

# Loop over selected cases
for index in case_indices:
    # Find current case name using index as identifies
    case_index = f"{index:03d}"
    case_name = [case for case in cases if case.startswith(f"case_{index:03d}")].pop()
    results_dir = os.path.join(SIMULATIONS_DIR, case_name, RESULTS_DIR)

    # Create the folder to save figures
    fig_dir = os.path.join(results_dir, "figures")
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    # Case parameters
    print(f"Processing {case_name}")
    labels = parse("case_{index}_{fluid}_T{temperature}K_P{pressure}bar", case_name)
    pressure = float(labels["pressure"].replace("_", "."))
    temperature = float(labels["temperature"].replace("_", "."))
    fluid = labels["fluid"]

    # Read residual file
    file_res = os.path.join(results_dir, f"case_{case_index}.trn")

    # Plot flow residuals
    fig, ax = plot_residuals(file_res)

    # Save plots
    filename = os.path.join(fig_dir, f"case_{case_index}_residuals")
    fig.savefig(filename + ".png", bbox_inches="tight")
    fig.savefig(filename + ".svg", bbox_inches="tight")
    fig.savefig(filename + ".pdf", bbox_inches="tight")

if show_figures:
    plt.show()
