import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import barotropy as bpy

bpy.set_plot_options(grid=True)


# Create folder to save results
DIR_OUT = "./output/validation"
os.makedirs(DIR_OUT, exist_ok=True)

# Location of the throat on the 3D CFD mesh
x_throat = -0.0259712
x_throat = -0.027

# Create the figure
fig, ax = plt.subplots(figsize=(6.0, 4.8))
ax.set_xlabel("Position along nozzle wall (mm)")
ax.set_ylabel("Static pressure (bar)")
ax.set_xlim([0, 80])
# ax.set_ylim([0, 70])

# Create the figure
fig1, ax1 = plt.subplots(figsize=(6.0, 4.8))
ax1.set_xlabel("Position along nozzle wall (mm)")
ax1.set_ylabel("Density (kg/m$^3$)")
ax1.set_xlim([0, 80])
# ax.set_ylim([0, 100])

# exp_case = "7b"
# exp_case = "4b"
exp_case = ["4b", "7b"]
# Load and plot experimental data
file_exp = "./nakagawa_2009_pressure_distributions.csv"
df_exp = pd.read_csv(file_exp)
df_Tdata = df_exp.loc[(df_exp["type"] == "Tdata") & (df_exp["case"].isin(exp_case))]
df_pdata = df_exp.loc[(df_exp["type"] == "pdata") & (df_exp["case"].isin(exp_case))]
df_IHEdata = df_exp.loc[(df_exp["type"] == "IHEdata") & (df_exp["case"].isin(exp_case))]
# df_Tdata = df_exp.loc[(df_exp["type"] == "Tdata") & (df_exp["case"] == exp_case)]
# df_pdata = df_exp.loc[(df_exp["type"] == "pdata") & (df_exp["case"] == exp_case)]
# df_IHEdata = df_exp.loc[(df_exp["type"] == "IHEdata") & (df_exp["case"] == exp_case)]
ax.plot(df_Tdata["x"]*1000, df_Tdata["p"]/1e5, "ko", markersize=4, label="From temperatures", zorder=3)
ax.plot(df_pdata["x"]*1000, df_pdata["p"]/1e5, "k^", markersize=4, label="From strain gauges",  zorder=3)
# ax.plot(df_IHEdata["x"]*1000, df_IHEdata["p"], color="darkorange", label="Isentropic HEM")

# Load and plot CFD data
DATAFILE = "./nakagawa_2009_cases.xlsx"
case_data = pd.read_excel(DATAFILE)
case_data = case_data[case_data["nozzle"] == "nozzle_2"]
# case_data = case_data[(case_data["tag"] == "equilibrium_sens_1") | (case_data["case"] == 500)]
case_data = case_data[(case_data["tag"] == "metastable_sens") | (case_data["tag"] == "roughness_sens")]

# cases = [1000, 1001, 1002, 1003]
# case_data = case_data[case_data["case"].isin(cases)]
case_data = case_data[np.isclose(case_data["q_onset"], 0.00) | (case_data["case"].isin([1000, 1001])) | (case_data["tag"] == "roughness_sens")]

# case_data = case_data[np.isclose(case_data["p0_in"], 61e5)]
# cases = [101, 102, 103, 104, 105]
# cases = [101, 102, 103, 104]
# case_data = case_data[case_data["case"].isin(cases)]
colors = plt.cm.magma(np.linspace(0.25, 0.75, len(case_data)))
for i, (index, row) in enumerate(case_data.iterrows()):
    case_name = f"case_{row["case"]}"
    dir_cfd = f"./output/{case_name}/latest/fluent_simulation"
    file_xy = os.path.join(dir_cfd, f"{case_name}_pressure.xy")
    df_cfd = bpy.fluent_automation.parse_fluent_xy(filename=file_xy)
    x_cfd = df_cfd["wall_centerline"]["Z-Coordinate"] - x_throat
    p_cfd = df_cfd["wall_centerline"]["Static Pressure"]

    file_xy = os.path.join(dir_cfd, f"{case_name}_density.xy")
    df_cfd = bpy.fluent_automation.parse_fluent_xy(filename=file_xy)
    d_cfd = df_cfd["wall_centerline"]["Density"]
    d_cfd = d_cfd[x_cfd < 0.077]
    p_cfd = p_cfd[x_cfd < 0.077]
    x_cfd = x_cfd[x_cfd < 0.077]
    # label = "CFD smooth wall" if row["case"] == 101 else rf"CFD $k_s=${row["wall_roughness_height"]*1e6:0.1f} $\mu$m"
    label = rf"CFD $Q_s=${row["q_onset"]:0.2f}, $Q_\text{{width}}=${row["q_transition"]:0.2f}"
    if row["case"] == "inviscid":
        ax.plot(x_cfd*1000, p_cfd/1e5, color="black", linestyle="--", label="CFD inviscid")
    else:
        ax.plot(x_cfd*1000, p_cfd/1e5, color=colors[i], label=label)

    if row["case"] == "inviscid":
        ax1.plot(x_cfd*1000, d_cfd, color="black", linestyle="--", label="CFD inviscid")
    else:
        ax1.plot(x_cfd*1000, d_cfd, color=colors[i], label=label)

ax.legend(loc="upper right", 
    handletextpad=0.4,  # Space between legend marker and text
    # columnspacing=0.5,  # Space between columns
    # borderaxespad=0.2,  # Padding between legend and plot area
    labelspacing=0.2,   # Vertical space between legend entries
    fontsize=9,
    ncols=1,
)
ax1.legend(loc="lower left", 
    handletextpad=0.4,  # Space between legend marker and text
    # columnspacing=0.5,  # Space between columns
    # borderaxespad=0.2,  # Padding between legend and plot area
    labelspacing=0.2,   # Vertical space between legend entries
    fontsize=10,
)

# Show figure
fig.tight_layout(pad=1)
bpy.savefig_in_formats(fig, os.path.join(DIR_OUT, "validation_3"))
plt.show()
