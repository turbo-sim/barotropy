# %% Import
import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI

# Sample data
FILE_IN = "berana_2009_pressure_profiles_original.xlsx"
FILE_OUT = "berana_2009_pressure_profiles_converted.xlsx"
data = pd.read_excel("berana_2009_data.xlsx")
pressure_data = pd.read_excel(FILE_IN)
fluidname = data["fluidname"].values[0]
pressure_data["pressure"] = pressure_data["pressure"] * 10**6


# Define conversion functions
def compute_converging_T(pressure, s_in):
    temperature = round(PropsSI("T", "P", pressure, "S", s_in, fluidname), 4)
    return temperature


def compute_diverging_T(pressure):
    temperature = round(PropsSI("T", "P", pressure, "Q", 0, fluidname), 4)
    return temperature


# Define a custom function that applies different calculations based on conditions
def convert_p2T(row):
    if row["type"] == "converging":
        # PropsSI(output, input1, value1, input2, value2, fluid)
        p_in = data.loc[data["Case"] == row["case"], "P_0_in"].values[0]
        T_in = data.loc[data["Case"] == row["case"], "T_0_in"].values[0] + 273.15
        s_in = PropsSI("S", "P", p_in, "T", T_in, fluidname)
        return compute_converging_T(row["pressure"], s_in)
    elif row["type"] == "diverging":
        return compute_diverging_T(row["pressure"])


# Apply the custom function row-wise to create column 'temperature'
pressure_data["temperature"] = pressure_data.apply(convert_p2T, axis=1)

with pd.ExcelWriter(FILE_OUT) as writer:
    pressure_data.to_excel(writer, sheet_name="Sheet1", index=False)

print(pressure_data)
