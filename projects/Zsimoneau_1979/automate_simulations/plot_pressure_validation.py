import os
from parse import parse
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from barotropy import parse_fluent_xy, set_plot_options


# Define plot settings
show_figures = True
set_plot_options(fontsize=14, grid=True)

# Define directory paths
SIMULATIONS_DIR = "../data_simulations/"
EXPERIMENTAL_DIR = "../data_experiments/"
RESULTS_DIR = "simulation_results"

# Define which cases to plot
# case_indices = [90, 91, 92, 93, 94, 95]
# case_indices = list(range(88, 96)) + list(range(100, 118))
# case_indices = list(range(88, 96)) + [100, 101]
case_indices = [92]

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

    # Create figure
    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.gca()
    ax.set_title(
        rf"Case {index}   |   $T_{{0,\mathrm{{in}}}}={temperature}$ K   |   $p_{{0,\mathrm{{in}}}}={pressure}$ bar"
    )
    ax.set_xlabel("Axial position (mm)")
    ax.set_ylabel("Static pressure (bar)")
    ax.set_xscale("linear")
    ax.set_yscale("linear")
    ax.set_ylim([0, np.ceil(pressure / 10) * 10])
    # ax.set_ylim([])
    # ax.set_xticks([])
    # ax.set_yticks([])

    # Plot simulation data
        # Load Fluent XY datafile
    file_sim = os.path.join(results_dir, f"case_{case_index}_pressure.xy")
    df_sim = parse_fluent_xy(file_sim)
    order = df_sim["wall"]["Position"].argsort()
    x = df_sim["wall"]["Position"][order] * 1e3
    y = df_sim["wall"]["Static Pressure"][order] / 1e5
    ax.plot(
        x,
        y,
        # color="black",
        linewidth=1.5,
        linestyle="-",
        marker="none",
        label="Simulation results",
    )

    # Plot experimental data
    file_exp = os.path.join(EXPERIMENTAL_DIR, case_name + ".csv")
    df_exp = pd.read_csv(file_exp, delimiter=",")
    ax.plot(
        df_exp["Position"] * 1e3,
        df_exp["Static Pressure"] / 1e5,
        marker="o",
        linewidth=1.25,
        # color="red",
        linestyle="None",
        label="Experimental data",
    )

    # Create legend
    leg = ax.legend(loc="upper right")

    # Adjust PAD
    fig.tight_layout(pad=1, w_pad=None, h_pad=None)

    # Save plots
    filename = os.path.join(fig_dir, f"case_{case_index}_validation")
    fig.savefig(filename + ".png", bbox_inches="tight")
    fig.savefig(filename + ".svg", bbox_inches="tight")
    fig.savefig(filename + ".pdf", bbox_inches="tight")





    # Create figure
    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.gca()
    # ax.set_title(
    #     rf"Case {index}   |   $T_{{0,\mathrm{{in}}}}={temperature}$ K   |   $p_{{0,\mathrm{{in}}}}={pressure}$ bar"
    # )
    ax.set_xlabel("Axial position (mm)")
    ax.set_ylabel("Density (kg/m$^3$)")
    ax.set_xscale("linear")
    ax.set_yscale("linear")
    file_sim = os.path.join(results_dir, f"case_{case_index}_density.xy")
    df_sim = parse_fluent_xy(file_sim)
    ax.plot(
        df_sim["wall"]["Position"],
        df_sim["wall"]["Density"],
        # color="black",
        linewidth=1.25,
        linestyle="-",
        marker="none",
        label="Fluent density",
    )
    file_sim = os.path.join(results_dir, f"case_{case_index}_barotropic_density.xy")
    df_sim_barotropic = parse_fluent_xy(file_sim)
    ax.plot(
        df_sim_barotropic["wall"]["Position"],
        df_sim_barotropic["wall"]["rhomass_barotropic"],
        # color="black",
        linewidth=1.25,
        linestyle="--",
        marker="none",
        label="Barotropic density",
    )
    ax.legend(loc="upper right")
    fig.tight_layout(pad=1, w_pad=None, h_pad=None)
    filename = os.path.join(fig_dir, f"case_{case_index}_density")
    fig.savefig(filename + ".png", bbox_inches="tight")
    fig.savefig(filename + ".svg", bbox_inches="tight")
    fig.savefig(filename + ".pdf", bbox_inches="tight")




    # Create figure
    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.gca()
    ax.set_xlabel("Axial position (mm)")
    ax.set_ylabel("Density deviation (%)")
    ax.set_xscale("linear")
    ax.set_yscale("linear")

    # Plot simulation data
    ax.plot(
        df_sim["wall"]["Position"],
        (df_sim["wall"]["Density"]/df_sim_barotropic["wall"]["rhomass_barotropic"]-1)*100,
        # color="black",
        linewidth=1.25,
        linestyle="-",
        marker="none",
        label="Fluent density",
    )
    filename = os.path.join(fig_dir, f"case_{case_index}_rhodev")
    fig.savefig(filename + ".png", bbox_inches="tight")
    fig.savefig(filename + ".svg", bbox_inches="tight")
    fig.savefig(filename + ".pdf", bbox_inches="tight")



if show_figures:
    plt.show()
