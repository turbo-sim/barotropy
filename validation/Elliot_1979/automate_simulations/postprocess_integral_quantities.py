import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd

from barotropy import parse_fluent_out

# Import data about simulation cases
case_data = pd.read_excel("simulation_cases.xlsx", skiprows=0)
# case_data = case_data[case_data["index"].isin([1])]
case_data = case_data[case_data["tag"].isin(["a", "b", "c", "d", "e"])]

# Define directory paths
SIMULATIONS_ROOTDIR = "output"
TARGET_DIR = "simulation_latest"
OUTPUT_FILE = "simulation_results.xlsx"

# Get values from the last row of the Fluent report files (.out)
last_rows = []

# Loop over selected cases
for _, row in case_data.iterrows():
    # Get the current case name from the Identifier and Case columns of Excel file
    index = int(row["index"])
    tag = row["tag"]
    case_name = f"case_{index:03d}" if pd.isna(tag) else f"case_{index:03d}_{tag}"
    results_dir = os.path.join(SIMULATIONS_ROOTDIR, case_name, TARGET_DIR)
    print(f"Processing {case_name}")

    # Load Fluent report file
    report_files = glob.glob(os.path.join(results_dir, "*.out"))
    if len(report_files) == 1:
        df_sim = parse_fluent_out(report_files[0])
    elif len(report_files) == 0:
        raise FileNotFoundError(f"No .out file found for {case_name} in {results_dir}")
    else:
        raise ValueError(f"Multiple .out files found for {case_name} in {results_dir}. Expected only one.")
    
    
    # Create a new series with Identifier, Case, and last_row
    data = {'index': index, 'tag': tag}
    data.update(df_sim.iloc[-1].to_dict())
    combined_series = pd.Series(data)

    # Append this combined series to the list
    last_rows.append(combined_series)

# Convert list of Series to dataframe
result_df = pd.DataFrame(last_rows).reset_index(drop=True)

# Export the results to excel
file_output = os.path.join(SIMULATIONS_ROOTDIR, OUTPUT_FILE)
result_df.to_excel(file_output, index=False)
print(result_df)




# fig = plt.figure(figsize=(6.4, 4.8))
# ax = fig.gca()
# ax.set_xscale("linear")
# ax.set_yscale("linear")

# # Import data
# ax.plot(
#     df_sim["Iteration"],
#     df_sim["velocity_outlet"],
#     # color="black",
#     linewidth=1.25,
#     linestyle="-",
#     marker="+",
#     label="Pressure outlet",
# )

# # ax.plot(
# #     df_sim["Iteration"],
# #     -df_sim["mass_flow_rate_throat"],
# #     # color="black",
# #     linewidth=1.25,
# #     linestyle="-",
# #     marker="none",
# #     label="Throat",
# # )

# # # Create legend
# # ax.legend(loc="upper right", fontsize=11)

# # # Adjust PAD
# # fig.tight_layout(pad=1, w_pad=None, h_pad=None)

# # plt.show()

