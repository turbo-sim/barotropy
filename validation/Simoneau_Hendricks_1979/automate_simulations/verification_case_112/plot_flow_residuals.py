import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd

from barotropy import set_plot_options, read_residual_file


# Define plot settings
show_figures = True
set_plot_options(fontsize=14, grid=True, margin=0.00)

# Create the folder to save figures
fig_dir = os.path.join("figures")
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

# Read residual file
file_res = os.path.join("./data_simulations/case_112_manual.trn")
df = read_residual_file(file_res, res_number=5)

# Create figure
fig = plt.figure(figsize=(7, 4))
ax = fig.gca()
ax.set_xlabel("Iteration number")
ax.set_ylabel("Normalized RMS residual")
ax.set_xscale("linear")
ax.set_yscale("log")
ax.set_xlim([0, 1e4])
# ax.set_ylim([])
# ax.set_xticks([])
# ax.set_yticks([])

# Plot the residuals
name_map = {
    "continuity": "$p$ - continuity",
    "x-velocity": "$v_x$ - momentum",
    "y-velocity": "$v_y$ - momentum",
    "k": "$k$ - turbulence",
    "omega": "$\omega$ - turbulence",
}

for column in df.columns[1:]:
    ax.plot(
        df["iter"],
        df[column],
        # color="black",
        linewidth=0.5,
        linestyle="-",
        marker="none",
        label=name_map[column],
    )

# Create legend
leg = ax.legend(loc="upper right", fontsize=7.5)

# Adjust PAD
fig.tight_layout(pad=1, w_pad=None, h_pad=None)

# Save plots
filename = os.path.join(fig_dir, f"case_112_residuals")
fig.savefig(filename + ".png", bbox_inches="tight")
fig.savefig(filename + ".svg", bbox_inches="tight")
fig.savefig(filename + ".pdf", bbox_inches="tight")
if show_figures:
    plt.show()
