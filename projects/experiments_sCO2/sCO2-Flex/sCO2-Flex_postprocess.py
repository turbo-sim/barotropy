import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd

import barotropy as bpy
import barotropy.fluids.core as props
import barotropy.sCO2_utilities as sco2


# ---------------------------------------------------------------------------- #
# Load and postprocess data
# ---------------------------------------------------------------------------- #

# Define case parameters
CASE = 'sCO2-Flex'
SAVE_FIGS = True
OUTPUT_DIR = "results"
bpy.set_plot_options(minor_ticks=False)

# Create output directory if it does not exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Import data about simulation cases
df = pd.read_excel(
    f"{CASE}_data_digitized.xlsx",
    sheet_name="data",
    skiprows=lambda x: x in [1],  # Skip unit row
)

# Extract the desired dataset
df = df[(df["compressor"] == CASE) & (df["data_source"] == "EXP")]

# Create fluid
fluid = props.Fluid(name="CO2", backend="HEOS", exceptions=True)

# Postprocess data
df = df.assign(**df.apply(lambda row: sco2.compute_inlet_state(row, fluid), axis=1))

# Save postprocessed data
df.to_excel(f"{CASE}_data_postprocessed.xlsx")


