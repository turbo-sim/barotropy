import os
from parse import parse
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd


from barotropy import parse_fluent_xy, set_plot_options, parse_fluent_out


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

# Get the case names
cases = os.listdir(SIMULATIONS_DIR)

# Import excel file with experimental mass flowrates
df_exp = pd.read_excel("../cases_summary.xlsx", skiprows=0)

# Initialize lists to store mass flow rates
mass_flow_sim = []
mass_flow_exp = []

# Loop over selected cases
for index in case_indices:
    # Find current case name using index as identifies
    case_index = f"{index:03d}"
    case_name = [case for case in cases if case.startswith(f"case_{index:03d}")].pop()
    results_dir = os.path.join(SIMULATIONS_DIR, case_name, RESULTS_DIR)

    # Load Fluent report file
    file_sim = os.path.join(results_dir, f"case_{case_index}.out")
    df_sim = parse_fluent_out(file_sim)

    # Get the simulated and experimental mass flow rates
    mass_flow_sim.append(df_sim["mass_in"].iloc[-1])
    mass_flow_exp.append(df_exp["m_exp"][df_exp["Case"] == index].iloc[0])


# Print comparison table
print(
    f"{'Case':>10} {'Mass flow exp. (kg/s)':>25} {'Mass flow sim. (kg/s)':>25} {'Abs. deviation (kg/s)':>25} {'Rel. deviation (%)':>25}"
)
for index, m_sim, m_exp in zip(case_indices, mass_flow_sim, mass_flow_exp):
    print(
        f"{index:>10} {m_sim:>25.6f} {m_exp:>25.6f} {m_sim-m_exp:>25.6f} {(m_sim-m_exp)/m_exp*100:>25.6f}"
    )

# Specify the output file path
output_file_path = os.path.join(SIMULATIONS_DIR, "mass_flow_comparison.md")

# Open the file for writing
with open(output_file_path, "w") as f:
    # Print the header to the file
    header = "| Case | Mass flow exp. (kg/s) | Mass flow sim. (kg/s) | Abs. deviation (kg/s) | Rel. deviation (%) |"
    f.write(header + "\n")
    f.write(
        "|------|-----------------------|-----------------------|-----------------------|--------------------|\n"
    )

    # Print data rows to the file
    for index, m_sim, m_exp in zip(case_indices, mass_flow_sim, mass_flow_exp):
        row = f"| {index:>4} | {m_exp:>20.6f} | {m_sim:>20.6f} | {(m_sim-m_exp):>20.6f} | {(m_sim-m_exp)/m_exp*100:>18.6f} |"
        f.write(row + "\n")

print(f"Table written to {output_file_path}")


# Make comparison plot
fig = plt.figure(figsize=(6.4, 5.8))
ax = fig.gca()
ax.set_title("Comparison of experimental and simulation data")
ax.set_xlabel("Mass flow rate from experiments (kg/s)")
ax.set_ylabel("Mass flow rate from simulations (kg/s)")
ax.set_xscale("linear")
ax.set_yscale("linear")
ax.set_xlim([0.0, 0.5])
ax.set_ylim([0.0, 0.5])

# Plot the lines and markers
x = np.linspace(0.0, 0.5, 100)
ax.plot(mass_flow_exp, mass_flow_sim, linestyle="none", marker="o")
ax.plot(
    x,
    1.00 * x,
    color="black",
    linestyle="-",
    label="$\dot{m}_\mathrm{sim} = \dot{m}_\mathrm{exp}$",
)
ax.plot(
    x,
    0.95 * x,
    color="black",
    linestyle="--",
    label="$\dot{m}_\mathrm{sim} = 0.95\, \dot{m}_\mathrm{exp}$",
)
ax.plot(
    x,
    1.05 * x,
    color="black",
    linestyle="--",
    label="$\dot{m}_\mathrm{sim} = 1.05\, \dot{m}_\mathrm{exp}$",
)

# Add a legend
plt.legend(loc="upper left")

# Adjust PAD
fig.tight_layout(pad=1, w_pad=None, h_pad=None)

# Save plots
filename = os.path.join(SIMULATIONS_DIR, "mass_flow_rate_comparison")
fig.savefig(filename + ".png", bbox_inches="tight")
fig.savefig(filename + ".svg", bbox_inches="tight")
fig.savefig(filename + ".pdf", bbox_inches="tight")

# Show the plot
if show_figures:
    plt.show()
