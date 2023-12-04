import os
import pandas as pd

# Define conversion functions
def kelvin_to_celsius(kelvin):
    return kelvin - 273.15

def pascal_to_bar(pascal):
    return pascal / 1e5

def unit_to_kilo(pascal):
    return pascal / 1e3

def no_conversion(x):
    return x

def psia_to_bar(psia):
    return psia * 0.0689476

def fahrenheit_to_celsius(fahrenheit):
    return (fahrenheit - 32) * 5 / 9


# Define the name of the common excel file
FILENAME = "sCO2_experimental_data.xlsx"

# Define which cases will be included in the common excel file
CASES = ["sCO2-HeRo", "sCO2-flex", "BMPC", "Sandia", "KAIST", "MIT"]

# Define which variables with be included in the common excel file
conversion_map = {
    "compressor": {'func': no_conversion, 'unit': '-'},
    "data_source": {'func': no_conversion, 'unit': '-'},
    "tag": {'func': no_conversion, 'unit': '-'},
    "T_in": {'func': kelvin_to_celsius, 'unit': 'degC'},
    "p_in": {'func': pascal_to_bar, 'unit': 'bar'},
    "rho_in": {'func': no_conversion, 'unit': 'kg/m3'},
    "Z_in": {'func': no_conversion, 'unit': '-'},
    "a_in": {'func': no_conversion, 'unit': 'm/s'},
    "s_in": {'func': unit_to_kilo, 'unit': 'kJ/kg/K'},
    "p_out": {'func': pascal_to_bar, 'unit': 'bar'},
    "pressure_ratio": {'func': no_conversion, 'unit': '-'},
    "isentropic_work": {'func': no_conversion, 'unit': 'kJ/kg'},
    "angular_speed": {'func': no_conversion, 'unit': 'RPM'},
    "mass_flow": {'func': no_conversion, 'unit': 'kg/s'},
    "efficiency": {'func': no_conversion, 'unit': '%'},
}

# Extract column names for reindexing
columns_to_export = list(conversion_map.keys())

# Initialize an empty list to store the dataframes
dfs = []

for case in CASES:
    # Read the excel file
    file = os.path.join(f"{case}", f"{case}_data_postprocessed.xlsx")
    df = pd.read_excel(file)

    # Select only the columns you want to export
    # Use reindex to add missing columns as NaNs
    df_selected = df.reindex(columns=conversion_map)

    # Apply unit conversions
    for column, info in conversion_map.items():
        if column in df_selected.columns:
            df_selected[column] = df_selected[column].apply(info['func'])

    # Append the selected columns to the list
    dfs.append(df_selected)

# Concatenate all the dataframes in the list
df_all = pd.concat(dfs, ignore_index=True)

# Insert a row with units at the beginning of the dataframe
units_row = {col: info['unit'] for col, info in conversion_map.items()}
df_all = pd.concat([pd.DataFrame([units_row]), df_all], ignore_index=True)

# Export the consolidated dataframe to Excel
df_all.to_excel(FILENAME, index=False)

