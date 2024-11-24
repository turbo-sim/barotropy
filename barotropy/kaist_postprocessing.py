import copy
import numpy as np
import pandas as pd

from . import fluid_properties as props


DATA_MAPPING = {
    "DPT-101 (bar)": "dP_recuperator_hot",
    "DPT-102 (bar)": "dP_cooler_PCHE",
    "DPT-103 (bar)": "dP_cooler_STHE",
    "DPT-104 (bar)": "dP_recuperator_cold",
    "DPT-301 (bar)": None,
    "DPT-302 (bar)": "dP_cooler_STHE_water",
    "FT-101D (kg/L)": "D_compressor_in",
    "FT-101F (kg/min)": "MF_compressor_in",
    "FT-102D (kg/L)": "D_compressor_out",
    "FT-102F (kg/min)": "MF_compressor_out",
    "FT-103F (kg/min)": "MF_turbine_in",
    "FT-301F (kg/min)": "MF_chiller_1",
    "FT-302F (kg/min)": "MF_chiller_2",
    "FT-303F (kg/min)": "MF_coolers",
    "FU,1 (hz)": "angular_speed",
    "FU,2 (hz)": None,
    "FU,3 (hz)": None,
    "FU,SIGM (hz)": None,
    "HE-001M (%)": None,
    "IRMS,1 (A)": None,
    "IRMS,2 (A)": None,
    "IRMS,3 (A)": None,
    "IRMS,SIGM (A)": None,
    "MFC-101M (%)": None,
    "MFC-101P (%)": None,
    "MFC-102M (%)": None,
    "MFC-102P (%)": None,
    "MFC-103M (%)": None,
    "MFC-103P (%)": None,
    "MFC-301M (%)": None,
    "MFC-301P (%)": None,
    "MFC-VENTM (%)": None,
    "MFC-VENTP (%)": None,
    "NONE (NONE)": None,
    "NONE (NONE) 1": None,
    "NONE (NONE) 2": None,
    "NONE (NONE) 3": None,
    "P,1 (W)": None,
    "P,2 (W)": None,
    "P,3 (W)": None,
    "P,SIGM (W)": None,
    "PT-101 (bar)": "P_turbine_in",
    "PT-102 (bar)": "P_turbine_out",
    "PT-103 (bar)": "P_recuperator_hot_in",
    "PT-105 (bar)": "P_recuperator_hot_out",
    "PT-107 (bar)": "P_cooler_PCHE_out",
    "PT-109 (bar)": "P_filter_in",
    "PT-109A (bar)": "P_compressor_in",
    "PT-110 (bar)": "P_compressor_out",
    "PT-111 (bar)": "P_recuperator_cold_in",
    "PT-113 (bar)": "P_heater_in",
    "PT-114 (bar)": "P_heater_out",
    "PT-201 (bar)": "P_cooler_PCHE_water_out",
    "PT-202 (bar)": "P_cooler_PCHE_water_in",
    "S,1 (VA)": None,
    "S,2 (VA)": None,
    "S,3 (VA)": None,
    "S,SIGM (VA)": None,
    "TE-101 (degC)": "T_turbine_in",
    "TE-102 (degC)": "T_turbine_out",
    "TE-103 (degC)": "T_recuperator_hot_in",
    "TE-104 (degC)": "T_recuperator_hot_out",
    "TE-105 (degC)": "T_cooler_PCHE_in",
    "TE-106 (degC)": "T_cooler_PCHE_out",
    "TE-107 (degC)": "T_cooler_STHE_in",
    "TE-108 (degC)": "T_cooler_STHE_out",
    "TE-109 (degC)": "T_filter_in",
    "TE-110 (degC)": "T_compressor_out",
    "TE-111 (degC)": "T_recuperator_cold_in",
    "TE-112 (degC)": "T_recuperator_cold_out",
    "TE-113 (degC)": "T_heater_in",
    "TE-114 (degC)": "T_heater_out",
    "TE-201 (degC)": None,
    "TE-202 (degC)": None,
    "TE-203 (degC)": "T_compressor_core",
    "TE-204 (degC)": "T_turbine_core",
    "TE-301 (degC)": "T_precooler_water_out",
    "TE-302 (degC)": "T_precooler_water_in",
    "TE-303 (degC)": "T_cooler_water_in",
    "TE-304 (degC)": "T_cooler_water_out",
    "URMS,1 (V)": None,
    "URMS,2 (V)": None,
    "URMS,3 (V)": None,
    "URMS,SIGM (V)": None,
}




def load_and_convert_data(data_file, data_mapping):
    # Load the data from the Excel file
    df = pd.read_excel(data_file)
    df_1 = pd.DataFrame(
        {col.replace(" mean", ""): df[col] for col in df.columns if " mean" in col}
    )
    df_2 = pd.DataFrame(
        {
            col.replace(" std_dev_mean", ""): df[col]
            for col in df.columns
            if " std_dev_mean" in col
        }
    )

    # Filter out None values from data_mapping dictionary
    filtered_mapping = {k: v for k, v in data_mapping.items() if v is not None}

    # Renaming the columns using the filtered mapping
    df_1.rename(columns=filtered_mapping, inplace=True)
    df_2.rename(columns=filtered_mapping, inplace=True)

    # Insert the 'tag' column as the first column, if it exists
    if "tag" in df.columns:
        df_1.insert(0, "tag", df["tag"])
        df_2.insert(0, "tag", df["tag"])

    # Convert units to SI system
    df_1 = convert_units(df_1, mode="absolute") 
    df_2 = convert_units(df_2, mode="uncertainty")

    return df_1, df_2


def convert_units(df, mode="absolute"):
    """
    Convert measured quantities to SI units.

    The function assumes:
    - Temperatures are in degrees Celsius and are converted to Kelvin.
    - Pressures are in gauge bars and are converted to Pascals (absolute).
    - Angular speeds are in Hz and are converted to rad/s.
    - Mass flows are in kg/min and are converted to kg/s.

    For 'absolute' mode, multiplicative and additive conversions are applied
    For 'uncertainty' mode, only multiplicative conversions are applied (typical deviation is a differential quantity)

    Parameters
    ----------
    df : DataFrame
        The pandas DataFrame containing the data to be converted.
    mode : str, optional
        Specifies whether the data in the DataFrame are 'absolute' measurements or 'uncertainty' values.
        Default is 'absolute'.

    Returns
    -------
    DataFrame
        The DataFrame with converted units.

    """

    # Work on a copy
    df = copy.copy(df)

    # Convert temperatures from Celsius to Kelvin (only if dealing with absolute values)
    if mode == "absolute":
        temp_cols = df.filter(regex="^T_").columns
        if not temp_cols.empty:
            df[temp_cols] = df[temp_cols] + 273.15

    # Convert pressures from gauge bar to Pascal
    pressure_cols = df.filter(regex="^P_").columns
    if not pressure_cols.empty:
        if mode == "absolute":
            df[pressure_cols] = df[pressure_cols] * 1e5 + 101325
        elif mode == "uncertainty":
            df[pressure_cols] = df[pressure_cols] * 1e5
        else:
            raise ValueError(
                f"Invalid value for 'mode': {mode}. Valid options are 'absolute' or 'uncertainty'"
            )

    # Convert angular speed from Hz to rad/s
    omega_cols = df.filter(regex="^angular_speed").columns
    if not omega_cols.empty:
        df[omega_cols] = df[omega_cols] * 2 * np.pi

    # Convert mass flow from kg/min to kg/s
    mass_flow_cols = df.filter(regex="^MF_").columns
    if not mass_flow_cols.empty:
        df[mass_flow_cols] = df[mass_flow_cols] / 60

    return df


def calculate_enthalpy_rise(
    T_filter_in, P_filter_in, P_filter_out, T_out, p_out, MF_in, MF_out, angular_speed
):
    # Create fluid object
    fluid_name = "CO2"
    fluid = props.Fluid(name=fluid_name, backend="HEOS", exceptions=True)

    # Assume saturated vapor if point is subcooled liquid due to measurement uncertainty
    if P_filter_in < fluid.critical_point.p:
        state_sat = fluid.get_state(props.PQ_INPUTS, P_filter_in, 1.00)
        T_filter_in = state_sat.T + 1e-3 if (T_filter_in < state_sat.T) else T_filter_in

    # Set states
    state_filter_in = fluid.get_state(props.PT_INPUTS, P_filter_in, T_filter_in)
    state_in = fluid.get_state(props.HmassP_INPUTS, state_filter_in.h, P_filter_out)
    state_out = fluid.get_state(props.PT_INPUTS, p_out, T_out)
    state_out_s = fluid.get_state(props.PSmass_INPUTS, p_out, state_in.s)

    # Calculate work and ideal work
    work = state_out.h - state_in.h  ## / state_in.a**2
    ideal_work = state_out_s.h - state_in.h
    efficiency = ideal_work / work

    # Reduced mass flow
    work_eq = (state_out.h - state_in.h) / state_filter_in.a**2
    MF_in_eq = MF_in / state_filter_in.d / state_filter_in.a

    # rho_in = state_in.rho
    # a_in = state_in.a
    # s_in = state_in.s
    # h_in = state_in.h
    # T_in = state_in.T
    # rho_in = state_in.rho
    a_in = state_filter_in.a
    s_in = state_filter_in.s
    h_in = state_filter_in.h
    T_in = state_filter_in.T
    p_in = state_filter_in.p
    rho_in = state_filter_in.rho

    return work, work_eq, MF_in_eq, h_in, T_in, p_in, rho_in, a_in, s_in
