import pandas as pd
import numpy as np

from . import properties as props

from typing import Dict, Tuple


PROPERTY_LIST = ["p", "T", "rho", "gamma", "Z", "a", "h", "s"]


def compute_inlet_state(row, fluid):

    # Check for existence of columns and non-NaN values, then calculate accordingly
    if "pressure_in" in row and "temperature_in" in row and not pd.isna(row["pressure_in"]) and not pd.isna(row["temperature_in"]):
        p_in = row["pressure_in"] * 1e5  # bar to Pa
        T_in = row["temperature_in"] + 273.15  # degC to K
        fluid.set_state(props.PT_INPUTS, p_in, T_in)
    elif "pressure_in" in row and "density_in" in row and not pd.isna(row["pressure_in"]) and not pd.isna(row["density_in"]):
        d_in = row["density_in"]
        p_in = row["pressure_in"] * 1e5  # bar to Pa
        fluid.set_state(props.DmassP_INPUTS, d_in, p_in)
    elif "density_in" in row and "temperature_in" in row and not pd.isna(row["density_in"]) and not pd.isna(row["temperature_in"]):
        d_in = row["density_in"]
        T_in = row["temperature_in"] + 273.15  # degC to K
        fluid.set_state(props.DmassT_INPUTS, d_in, T_in)
    elif "pressure_in" in row and "entropy_in" in row and not pd.isna(row["pressure_in"]) and not pd.isna(row["entropy_in"]):
        s_in = row["entropy_in"] * 1e3  # kJ/kg/K to J/kg/K
        p_in = row["pressure_in"] * 1e5  # bar to Pa
        fluid.set_state(props.PSmass_INPUTS, p_in, s_in)
    elif "temperature_in" in row and "entropy_in" in row and not pd.isna(row["temperature_in"]) and not pd.isna(row["entropy_in"]):
        s_in = row["entropy_in"] * 1e3  # kJ/kg/K to J/kg/K
        T_in = row["temperature_in"] + 273.15  # degC to K
        fluid.set_state(props.SmassT_INPUTS, s_in, T_in)
    else:
        raise ValueError("Inlet state must be specified as 'pressure-temperature', 'pressure-density', or 'temperature-density' and the respective columns must exist.")

    # Append '_in' to property names and return
    inlet_props = {
        key + "_in": value
        for key, value in fluid.properties.items()
        if key in PROPERTY_LIST
    }

    return pd.Series(inlet_props)


def compute_outlet_state_isentropic(row, fluid):
    # Calculate inlet properties with p-T function call
    p_in = row["pressure_in"] * 1e5  # bar to Pa
    T_in = row["temperature_in"] + 273.15  # degC to K
    state_in = fluid.set_state(props.PT_INPUTS, p_in, T_in)

    # Calculate isentropic exit state
    p_out = row["pressure_ratio"] * p_in
    state_out = fluid.set_state(props.PSmass_INPUTS, p_out, state_in.s)

    # Append '_out' to exit properties and include dh_is
    outlet_props = {
        key + "_out": value
        for key, value in state_out.to_dict().items()
        if key in PROPERTY_LIST
    }
    outlet_props["isentropic_work"] = (state_out.h - state_in.h)/1e3

    return pd.Series(outlet_props)


def compute_critical_state(row, fluid):
    # Calculate inlet properties with p-T function call
    p_in = row["pressure_in"] * 1e5  # bar to Pa
    T_in = row["temperature_in"] + 273.15  # degC to K

    # Calculate critical state (i.e., sonic state for an isentropic expansion)
    state_critical = fluid.compute_sonic_state(props.PT_INPUTS, p_in, T_in)

    # Append '_out' to property names
    critical_props = {
        key + "_crit": value
        for key, value in state_critical.to_dict().items()
        if key in PROPERTY_LIST
    }

    return pd.Series(critical_props)


def calculate_scaling_variables(
    props: Dict[str, float]
) -> Tuple[float, float, float, float, float]:
    p = props["p"]
    T = props["T"]
    Z = props["Z"]  # Setting 1 here yields bad matching for HeRo compressor
    g = props["gamma"]
    v = np.sqrt(2 * g / (g + 1) * Z * T)
    f = g * (2 / (g + 1)) ** (g / (g - 1.0))
    a = props["a"]
    rho = props["rho"]
    return p, v, f, a, rho


def compute_equivalent_from_actual(
    row: pd.Series,
    state_ref: Dict[str, float],
    method: str = "perfect_gas",
) -> pd.Series:
    """Compute equivalent variables from actual conditions and reference state"""

    # Common calculations
    state_in = {col.rstrip("_in"): row[col] for col in row.index if col.endswith("_in")}
    p_1, v_1, f_1, a_1, rho_1 = calculate_scaling_variables(state_in)
    p_2, v_2, f_2, a_2, rho_2 = calculate_scaling_variables(state_ref)

    # Compute equivalent variables based on the selected method
    variables = {}
    if method == "perfect_gas":
        delta = p_1 / p_2
        theta = (v_1 / v_2) ** 2
        epsilon = f_2 / f_1
        variables["delta"] = delta
        variables["theta"] = theta
        variables["epsilon"] = epsilon
        variables["mass_flow_eq"] = row["mass_flow"] * (np.sqrt(theta) * epsilon / delta)
        variables["isentropic_work_eq"] = row["isentropic_work"] * (1 / theta)
        variables["angular_speed_eq"] = row["angular_speed"] * (1 / np.sqrt(theta))
    elif method == "real_gas":
        variables["mass_flow_eq"] = row["mass_flow"] * (rho_2 * a_2) / (rho_1 * a_1)
        variables["isentropic_work_eq"] = row["isentropic_work"] * (a_2 / a_1) ** 2
        variables["angular_speed_eq"] = row["angular_speed"] * (a_2 / a_1)

    return pd.Series(variables)


def compute_actual_from_equivalent(
    row: pd.Series,
    state_ref: Dict[str, float],
    method: str = "perfect_gas",
) -> pd.Series:
    """Compute equivalent variables from actual conditions and reference state"""

    # Common calculations
    state_in = {col.rstrip("_in"): row[col] for col in row.index if col.endswith("_in")}
    p_1, v_1, f_1, a_1, rho_1 = calculate_scaling_variables(state_in)
    p_2, v_2, f_2, a_2, rho_2 = calculate_scaling_variables(state_ref)

    # Compute equivalent variables based on the selected method
    variables = {}
    if method == "perfect_gas":
        delta = p_1 / p_2
        theta = (v_1 / v_2) ** 2
        epsilon = f_2 / f_1
        variables["delta"] = delta
        variables["theta"] = theta
        variables["epsilon"] = epsilon
        variables["mass_flow"] = row["mass_flow_eq"] * (np.sqrt(theta) * epsilon / delta) ** -1
        # variables["isentropic_work"] = row["isentropic_work_eq"] * (1 / theta) ** -1
        # variables["angular_speed"] = row["angular_speed_eq"] * (1 / np.sqrt(theta)) ** -1
    elif method == "real_gas":
        variables["mass_flow"] = row["mass_flow_eq"] * ((rho_2 * a_2) / (rho_1 * a_1))  ** -1
        # variables["isentropic_work"] = row["isentropic_work_eq"] * (a_2 / a_1) ** -2
        # variables["angular_speed"] = row["angular_speed_eq"] * (a_2 / a_1) ** -1

    return pd.Series(variables)


import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def plot_ps_diagram(fluid, ax, range_x, range_y, range_z):

    # Define nice colormap
    colormap = mcolors.LinearSegmentedColormap.from_list(
        "truncated_blues", plt.cm.Blues(np.linspace(0.4, 1.0, 256))
    )

    # Define limits
    prop_x = "s"
    prop_y = "p"
    prop_z = "T"
    prop_pair = props.PSmass_INPUTS

    # Plot phase diagram
    ax = fluid.plot_phase_diagram(
        prop_x,
        prop_y,
        plot_critical_point=True,
        plot_quality_isolines=True,
        plot_pseudocritical_line=True,
    )

    # Plot pressure isobars
    prop_dict = props.compute_property_grid(fluid, prop_pair, range_y, range_x)
    contour = ax.contour(
        prop_dict[prop_x],
        prop_dict[prop_y],
        prop_dict[prop_z],
        range_z,
        cmap=colormap,
        linewidths=0.5,
    )

    # Colorbar settings
    sm = plt.cm.ScalarMappable(
        cmap=colormap,
        norm=plt.Normalize(vmin=min(range_z), vmax=max(range_z)),
    )
    sm.set_array([])  # This line is necessary for ScalarMappable
    cbar = plt.colorbar(sm, ax=ax)
    cbar.ax.set_ylabel("Temperature [$^\circ$C]", rotation=90, labelpad=20)
    cbar.set_ticks(range_z)  # Set tick positions
    cbar.set_ticklabels([f"{level-273.15:.0f}" for level in range_z])

    return ax


def plot_Ts_diagram(fluid, ax, range_x, range_y, range_z):

    # Define nice colormap
    colormap = mcolors.LinearSegmentedColormap.from_list(
        "truncated_blues", plt.cm.Blues(np.linspace(0.4, 1.0, 256))
    )

    # Define limits
    prop_x = "s"
    prop_y = "T"
    prop_z = "p"
    prop_pair = props.SmassT_INPUTS
    
    # Plot phase diagram
    ax = fluid.plot_phase_diagram(
        prop_x,
        prop_y,
        plot_critical_point=True,
        plot_quality_isolines=True,
        plot_pseudocritical_line=True,
    )

    # Plot pressure isobars
    prop_dict = props.compute_property_grid(fluid, prop_pair, range_x, range_y)
    contour = ax.contour(
        prop_dict[prop_x],
        prop_dict[prop_y],
        prop_dict[prop_z],
        range_z,
        cmap=colormap,
        linewidths=0.5,
    )

    # Colorbar settings
    sm = plt.cm.ScalarMappable(
        cmap=colormap,
        norm=plt.Normalize(vmin=min(range_z), vmax=max(range_z)),
    )
    sm.set_array([])  # This line is necessary for ScalarMappable
    cbar = plt.colorbar(sm, ax=ax)
    cbar.ax.set_ylabel("Pressure [bar]", rotation=90, labelpad=20)
    cbar.set_ticks(range_z)  # Set tick positions
    cbar.set_ticklabels([f"{level/1e5:.0f}" for level in range_z])

    return ax
