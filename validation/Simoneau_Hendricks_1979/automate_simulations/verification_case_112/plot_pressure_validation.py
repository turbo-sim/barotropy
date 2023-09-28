import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from barotropy import parse_fluent_xy, set_plot_options


# Define plot settings
show_figures = True
set_plot_options(fontsize=14, grid=True)

# Create the folder to save figures
fig_dir = os.path.join("figures")
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

# Create figure
fig = plt.figure(figsize=(6.4, 4.8))
ax = fig.gca()
ax.set_xlabel("Axial position (mm)")
ax.set_ylabel("Static pressure (bar)")
ax.set_xscale("linear")
ax.set_yscale("linear")
ax.set_ylim([0, 60])
# ax.set_ylim([])
# ax.set_xticks([])
# ax.set_yticks([])


# Plot simulation data
sim_files = {
    "2000 iterations": "./data_simulations/case_112_pressure_2000iter.xy",
    "6500 iteration": "./data_simulations/case_112_pressure_6500iter.xy",
    "10000 iterations": "./data_simulations/case_112_pressure_8500iter.xy",
}

for label, file in sim_files.items():
    df_sim = parse_fluent_xy(file)
    order = df_sim["wall"]["Position"].argsort()
    x = df_sim["wall"]["Position"][order] * 1e3
    y = df_sim["wall"]["Static Pressure"][order] / 1e5
    ax.plot(
        x,
        y,
        # color="black",
        linewidth=1.0,
        linestyle="-",
        marker="none",
        label=label,
    )

# Plot experimental data
file_exp = os.path.join("./data_experiments/case_112_nitrogen_T129_70K_P56_20bar.csv")
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
filename = os.path.join(fig_dir, f"case_112_validation")
fig.savefig(filename + ".png", bbox_inches="tight")
fig.savefig(filename + ".svg", bbox_inches="tight")
fig.savefig(filename + ".pdf", bbox_inches="tight")

if show_figures:
    plt.show()
