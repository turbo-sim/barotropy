import os
import matplotlib.pyplot as plt

from .fluent_automation import parse_fluent_xy
from ..utilities import savefig_in_formats


def plot_nozzle_xy_data(case_name, dir_xy, data_vars, zone_list, x_var):
    
    # Parse all the requested files at the start
    parsed_data = {}
    for var in data_vars:
        filepath = os.path.join(dir_xy, f"{case_name}_{var}.xy")
        parsed_data[var] = parse_fluent_xy(filepath)

    # Define dictionary of variable/functions
    var_funcs = {
        "pressure": plot_pressure_distribution,
        "velocity": plot_velocity_distribution,
        "density": plot_density_distribution,
        "speed_sound": plot_speed_sound_distribution,
        "viscosity": plot_viscosity_distribution,
        "void_fraction": plot_void_fraction_distribution,
        "vapor_quality": plot_vapor_quality_distribution,
        "y_plus": plot_yplus_distribution,
    }

    # Loop only over the functions whose keys are substrings in any of the keys of data_vars
    figs = []
    for var_func_key in var_funcs:
        if any(var_func_key in data_var for data_var in data_vars):
            base_name = f"{case_name}_{var_func_key}"
            fig, _ = var_funcs[var_func_key](parsed_data, zone_list, x_var)
            savefig_in_formats(fig, os.path.join(dir_xy, base_name))
            figs.append(fig)
            
    return figs


def plot_pressure_distribution(parsed_data, zone_list, x_var):
    df = parsed_data.get("pressure")
    lines = []
    for zone in zone_list:
        add_line(lines, df, zone, "Static Pressure")

    return generic_plot_distribution(
        lines,
        title="Pressure distribution", 
        xlabel="Axial position (mm)", 
        ylabel="Static pressure (bar)", 
        scale_x=1e3,
        scale_y=1e-5,
        x_var=x_var,
    )

def plot_velocity_distribution(parsed_data, zone_list, x_var):
    df = parsed_data.get("velocity")
    lines = []
    for zone in zone_list:
        add_line(lines, df, zone, "Velocity Magnitude")

    return generic_plot_distribution(
        lines,
        title="Velocity distribution", 
        xlabel="Axial position (mm)", 
        ylabel="Fluid velocity (m/s)", 
        scale_x=1e3,
        scale_y=1,
        x_var=x_var,
    )

def plot_density_distribution(parsed_data, zone_list, x_var):
    df1 = parsed_data.get("density")
    df2 = parsed_data.get("barotropic_density")
    lines = []
    for zone in zone_list:
        add_line(lines, df1, zone, "Density", label="Fluent data")
        add_line(lines, df2, zone, "barotropic_density", linestyle="--", label="Barotropic model")

    return generic_plot_distribution(
        lines,
        title="Density distribution", 
        xlabel="Axial position (mm)", 
        ylabel="Density (kg/m$^3$)", 
        scale_x=1e3,
        scale_y=1,
        x_var=x_var,
    )

def plot_viscosity_distribution(parsed_data, zone_list, x_var):
    df1 = parsed_data.get("viscosity")
    df2 = parsed_data.get("barotropic_viscosity")
    lines = []
    for zone in zone_list:
        add_line(lines, df1, zone, "Molecular Viscosity", label="Fluent data")
        add_line(lines, df2, zone, "barotropic_viscosity", linestyle="--", label="Barotropic model")
    
    return generic_plot_distribution(
        lines,
        title="Viscosity distribution", 
        xlabel="Axial position (mm)", 
        ylabel=r"Viscosity ($\mathrm{\mu}$Pa s)", 
        scale_x=1e3,
        scale_y=1e+6,
        x_var=x_var,
    )

def plot_speed_sound_distribution(parsed_data, zone_list, x_var):
    df1 = parsed_data.get("speed_sound")
    df2 = parsed_data.get("barotropic_speed_sound")
    lines = []
    for zone in zone_list:
        add_line(lines, df1, zone, "Sound Speed", label="Fluent data")
        add_line(lines, df2, zone, "barotropic_speed_sound", linestyle="--", label="Barotropic model")

    return generic_plot_distribution(
        lines,
        title="Speed of sound distribution", 
        xlabel="Axial position (mm)", 
        ylabel="Speed of sound (m/s)", 
        scale_x=1e3,
        scale_y=1,
        x_var=x_var,
    )

def plot_void_fraction_distribution(parsed_data, zone_list, x_var):
    df = parsed_data.get("barotropic_void_fraction")
    lines = []
    for zone in zone_list:
        add_line(lines, df, zone, "barotropic_void_fraction", label="Barotropic model")

    return generic_plot_distribution(
        lines,
        title="Void fraction distribution", 
        xlabel="Axial position (mm)", 
        ylabel="Void fraction (%)", 
        scale_x=1e3,
        scale_y=100,
        x_var=x_var,
    )

def plot_vapor_quality_distribution(parsed_data, zone_list, x_var):
    df = parsed_data.get("barotropic_vapor_quality")
    lines = []
    for zone in zone_list:
        add_line(lines, df, zone, "barotropic_vapor_quality", label="Barotropic model")

    return generic_plot_distribution(
        lines,
        title="Mass fraction distribution", 
        xlabel="Axial position (mm)", 
        ylabel="Mass fraction (%)", 
        scale_x=1e3,
        scale_y=100,
        x_var=x_var,
    )

def plot_yplus_distribution(parsed_data, zone_list, x_var):
    df = parsed_data.get("y_plus")
    lines = []
    for zone in zone_list:
        add_line(lines, df, zone, "Wall Yplus")

    return generic_plot_distribution(
        lines,
        title="$y+$ distribution", 
        xlabel="Axial position (mm)", 
        ylabel="$y+$ at the wall", 
        scale_x=1e3,
        scale_y=1,
        x_var=x_var,
    )



# Helper functions
def add_line(lines, df, zone, data_key, label="Fluent data", linestyle="-"):
    if df is not None and data_key in df[zone].columns:
        lines.append({
            "df": df, 
            "zone": zone, 
            "data_key": data_key, 
            "legend": label, 
            "linestyle": linestyle
        })


def generic_plot_distribution(lines, fig=None, ax=None, title="", xlabel="", ylabel="", scale_x=1, scale_y=1, x_var="X-Coordinate"):
    if not fig and not ax:
        fig, ax = plt.subplots(figsize=(6.4, 4.8))
    
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    for line in lines:
        df = line.get('df')
        zone = line.get('zone', "")
        data_key = line.get('data_key', "")
        legend = line.get('legend', None)
        linestyle = line.get('linestyle', "-")
        
        ax.plot(df[zone][x_var] * scale_x, df[zone][data_key] * scale_y, linewidth=1.25, linestyle=linestyle, label=legend)
    
    if len(lines) > 1:
        ax.legend()
    
    fig.tight_layout(pad=1)

    return fig, ax
