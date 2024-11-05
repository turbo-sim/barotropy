import os
import matplotlib.pyplot as plt

from .fluent_automation import parse_fluent_xy
from ..utilities import savefig_in_formats


def plot_nozzle_data(case_name, fig_dir, data_vars):
    root_dir = os.path.dirname(fig_dir)
    
    # Parse all the requested files at the start
    parsed_data = {}
    for var in data_vars:
        filepath = os.path.join(root_dir, f"{case_name}_{var}.xy")
        parsed_data[var] = parse_fluent_xy(filepath)

    # Define dictionary of variable/functions
    var_funcs = {
        "pressure": plot_pressure_distribution,
        "velocity": plot_velocity_distribution,
        "density": plot_density_distribution,
        "speed_sound": plot_speed_sound_distribution,
        "viscosity": plot_viscosity_distribution,
        "void_fraction": plot_void_fraction_distribution,
        "mass_fraction": plot_mass_fraction_distribution,
        "y_plus": plot_yplus_distribution,
    }

    # Loop only over the functions whose keys are substrings in any of the keys of data_vars
    for var_func_key in var_funcs:
        if any(var_func_key in data_var for data_var in data_vars):
            base_name = f"{case_name}_{var_func_key}"
            fig, _ = var_funcs[var_func_key](parsed_data)
            savefig_in_formats(fig, os.path.join(fig_dir, base_name))


def plot_pressure_distribution(parsed_data):
    df = parsed_data.get("pressure")
    lines = []
    add_line(lines, df, "wall", "Static Pressure")
    return generic_plot_distribution(
        lines,
        title="Pressure distribution", 
        xlabel="Axial position (mm)", 
        ylabel="Static pressure (bar)", 
        scale_x=1e3,
        scale_y=1e-5
    )

def plot_velocity_distribution(parsed_data):
    df = parsed_data.get("velocity")
    lines = []
    add_line(lines, df, "axis", "Velocity Magnitude")
    return generic_plot_distribution(
        lines,
        title="Velocity distribution", 
        xlabel="Axial position (mm)", 
        ylabel="Fluid velocity (m/s)", 
        scale_x=1e3,
        scale_y=1
    )

def plot_density_distribution(parsed_data):
    df1 = parsed_data.get("density")
    df2 = parsed_data.get("barotropic_density")
    lines = []
    add_line(lines, df1, "wall", "Density", "Fluent data")
    add_line(lines, df2, "wall", "barotropic_density", "Barotropic model", linestyle="--")
    return generic_plot_distribution(
        lines,
        title="Density distribution", 
        xlabel="Axial position (mm)", 
        ylabel="Density (kg/m$^3$)", 
        scale_x=1e3,
        scale_y=1
    )

def plot_viscosity_distribution(parsed_data):
    df1 = parsed_data.get("viscosity")
    df2 = parsed_data.get("barotropic_viscosity")
    lines = []
    add_line(lines, df1, "wall", "Molecular Viscosity")
    add_line(lines, df2, "wall", "barotropic_viscosity", "Barotropic model", linestyle="--")
    return generic_plot_distribution(
        lines,
        title="Viscosity distribution", 
        xlabel="Axial position (mm)", 
        ylabel=r"Viscosity ($\mathrm{\mu}$Pa s)", 
        scale_x=1e3,
        scale_y=1e-6
    )

def plot_speed_sound_distribution(parsed_data):
    df1 = parsed_data.get("speed_sound")
    df2 = parsed_data.get("barotropic_speed_sound")
    lines = []
    add_line(lines, df1, "wall", "Sound Speed")
    add_line(lines, df2, "wall", "barotropic_speed_sound", "Barotropic model", linestyle="--")
    return generic_plot_distribution(
        lines,
        title="Speed of sound distribution", 
        xlabel="Axial position (mm)", 
        ylabel="Speed of sound (m/s)", 
        scale_x=1e3,
        scale_y=1
    )

def plot_void_fraction_distribution(parsed_data):
    df = parsed_data.get("barotropic_void_fraction")
    lines = []
    add_line(lines, df, "wall", "barotropic_void_fraction", "Barotropic model")
    return generic_plot_distribution(
        lines,
        title="Void fraction distribution", 
        xlabel="Axial position (mm)", 
        ylabel="Void fraction (%)", 
        scale_x=1e3,
        scale_y=100
    )

def plot_mass_fraction_distribution(parsed_data):
    df = parsed_data.get("barotropic_mass_fraction")
    lines = []
    add_line(lines, df, "wall", "barotropic_mass_fraction", "Barotropic model")
    return generic_plot_distribution(
        lines,
        title="Mass fraction distribution", 
        xlabel="Axial position (mm)", 
        ylabel="Mass fraction (%)", 
        scale_x=1e3,
        scale_y=100
    )

def plot_yplus_distribution(parsed_data):
    df = parsed_data.get("y_plus")
    lines = []
    add_line(lines, df, "wall", "Wall Yplus")
    return generic_plot_distribution(
        lines,
        title="$y+$ distribution", 
        xlabel="Axial position (mm)", 
        ylabel="$y+$ at the wall", 
        scale_x=1e3,
        scale_y=1
    )



# Helper functions
def add_line(lines, df, zone, data_key, legend="Fluent data", linestyle="-"):
    if df is not None and data_key in df[zone].columns:
        lines.append({
            "df": df, 
            "zone": zone, 
            "data_key": data_key, 
            "legend": legend, 
            "linestyle": linestyle
        })


def generic_plot_distribution(lines, fig=None, ax=None, title="", xlabel="", ylabel="", scale_x=1, scale_y=1):
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
        
        ax.plot(df[zone]["Position"] * scale_x, df[zone][data_key] * scale_y, linewidth=1.25, linestyle=linestyle, label=legend)
    
    if len(lines) > 1:
        ax.legend()
    
    fig.tight_layout(pad=1)

    return fig, ax
