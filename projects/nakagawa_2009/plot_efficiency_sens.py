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

# Create the figure
fig, ax = plt.subplots(figsize=(6.0, 4.8))
ax.set_xlabel("Position along nozzle wall (mm)")
ax.set_ylabel("Velocity (m/s)")
ax.set_xlim([0, 80])
# ax.set_ylim([0, 100e5])

# Load and plot CFD data
DATAFILE = "./nakagawa_2009_cases.xlsx"
case_data = pd.read_excel(DATAFILE)
case_data = case_data[case_data["nozzle"] == "nozzle_2"]
case_data = case_data[(case_data["tag"] == "efficiency_sens") | (case_data["tag"] == "inviscid")]
case_data = case_data[np.isclose(case_data["p0_in"], 91e5)]
# cases = [101, 102, 103, 104, 105]
# cases = [101, 102, 103, 104]
# case_data = case_data[case_data["case"].isin(cases)]
colors = plt.cm.magma(np.linspace(0.25, 0.75, len(case_data)))
for i, (index, row) in enumerate(case_data.iterrows()):
    case_name = f"case_{row["case"]}"
    print(case_name)
    dir_cfd = f"./output/{case_name}/latest/fluent_simulation"
    file_xy = os.path.join(dir_cfd, f"{case_name}_velocity.xy")
    df_cfd = bpy.fluent_automation.parse_fluent_xy(filename=file_xy)
    x_cfd = df_cfd["wall_centerline"]["Z-Coordinate"] - x_throat
    v_cfd = df_cfd["wall_centerline"]["Velocity Magnitude"]

    v_cfd = v_cfd[x_cfd < 0.077]
    x_cfd = x_cfd[x_cfd < 0.077]
    
    label = "CFD inviscid" if row["tag"] == "inviscid" else rf"CFD $\eta_p=${row["polytropic_efficiency"]*1e2:0.1f} %"
    ax.plot(x_cfd*1000, v_cfd, color=colors[i], label=label)


velocity_out = {}
for i, (index, row) in enumerate(case_data.iterrows()):
    case_name = f"case_{row["case"]}"
    dir_cfd = f"./output/{case_name}/latest/fluent_simulation"
    # print(case_name)
    file_out = os.path.join(dir_cfd, f"{case_name}.out")
    df = bpy.parse_fluent_out(file_out)

    velocity_out[case_name] = df["velocity_out"].iloc[-1]

print(velocity_out)

efficiency = (velocity_out["case_201"] / velocity_out["case_200"])**2

print(f"Nozzle efficiency: {100 * efficiency:0.2f} %")

# for i, (index, row) in enumerate(case_data.iterrows()):
#     case_name = f"case_{row["case"]}"
#     print(case_name)
#     dir_cfd = f"./output/{case_name}/latest/fluent_simulation"
#     file_xy = os.path.join(dir_cfd, f"{case_name}_velocity.xy")
#     df_cfd = bpy.fluent_automation.parse_fluent_xy(filename=file_xy)

#     try:
#         x_cfd = df_cfd["nozzle_centerline"]["Z-Coordinate"] - x_throat
#         v_cfd = df_cfd["nozzle_centerline"]["Velocity Magnitude"]

#         v_cfd = v_cfd[x_cfd < 0.077]
#         x_cfd = x_cfd[x_cfd < 0.077]
        
#         label = "CFD inviscid" if row["tag"] == "inviscid" else rf"CFD $\eta_p=${row["polytropic_efficiency"]*1e2:0.1f} %"
#         ax.plot(x_cfd*1000, v_cfd, "g--", label=label)
#     except:
#         pass



# ax.legend(loc="upper left", 
#     handletextpad=0.4,  # Space between legend marker and text
#     # columnspacing=0.5,  # Space between columns
#     # borderaxespad=0.2,  # Padding between legend and plot area
#     labelspacing=0.3,   # Vertical space between legend entries
#     fontsize=11,
# )

# # Show figure
# fig.tight_layout(pad=1)
# bpy.savefig_in_formats(fig, os.path.join(DIR_OUT, "efficiency_sens"))
# plt.show()
