import re
import os
import sys
import time
import shutil
import subprocess
import pexpect.popen_spawn

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from parse import parse

# TODO: This might be broken, verify that it works after I refactored into package
# from barotropy import stream_residuals, stream_transcript  # Import as module to get its path
from . import stream_residuals, stream_transcript  # Import as module to get its path
from ..utilities import is_float, wait_for_file


def create_fluent_journal(case_name, output_dir, template_filename, template_map, barotropic_vars=None, export_vars=None, export_zones=None, solution_strategy=None):
    """
    Create and validate a Fluent journal based on the provided parameters.

    This function generates a Fluent journal, which is a script containing a sequence 
    of commands for the Fluent simulation software. It uses a template file and a map 
    of variable values to generate this script, with options for adding additional 
    commands related to barotropic variables, export variables, export zones, and 
    a solution strategy.

    Parameters
    ----------
    case_name : str
        Name of the simulation case.
    output_dir : str
        Path to the directory of the current simulation case
    template_filename : str
        Path to the template file used to generate the journal.
    template_map : dict
        Mapping of variable names to their values used for rendering the template.
    barotropic_vars : list of str, optional
        List of barotropic variable expressions to be added to the journal.
    export_vars : dict, optional
        Dictionary mapping output filenames to Fluent variable names to be exported.
    export_zones : list of str, optional
        List of zones from which data should be exported.
    solution_strategy : dict of list, optional
        Dictionary where each value is a list of commands or options for the solution strategy.

    Returns
    -------
    str
        The generated Fluent journal as a string.
    """

    # Validate parameters
    if not output_dir:
        raise ValueError("Output directory must be provided.")
    elif not os.path.exists(output_dir):
        raise FileNotFoundError(f"Output directory {output_dir} does not exist.")

    if not template_filename:
        raise ValueError("Template filename must be provided.")
    elif not os.path.exists(template_filename):
        raise FileNotFoundError(f"Template file {template_filename} not found.")

    if export_zones and not export_vars:
        raise ValueError("export_vars must be provided if export_zones is specified.")
    
    if export_vars and (not export_zones or len(export_zones) == 0):
        raise ValueError("Export zones must not be empty if export_vars is provided.")

    if export_vars and not isinstance(export_vars, dict):
        raise ValueError("export_vars should be a dictionary mapping between export filenames and Fluent variable names.")

    if export_zones and not isinstance(export_zones, list):
        raise ValueError("export_zones should be a list of zones to export data from.")

    if solution_strategy and not (isinstance(solution_strategy, dict) and all(isinstance(v, list) for v in solution_strategy.values())):
        raise ValueError("solution_strategy should be a dictionary where each value is a list of commands or options.")
    
    # Path handling
    case_path = os.path.dirname(output_dir)
    file_prefix = os.path.join(output_dir, case_name)

    # Generate the journal file (a sequence of commands) for Fluent based on the template
    journal = render_fluent_journal(template_filename, template_map)

    # Add expressions for the barotropic model to the journal
    if barotropic_vars:
        path = os.path.join(case_path, "barotropic_model", "fluent_expressions.txt")
        journal = add_barotropic_expressions(journal, path, barotropic_vars)

    # Add solution strategy commands to the journal
    if solution_strategy:
        journal = add_solution_strategy(journal, solution_strategy)

    # Add commands to export specific variables after the simulation
    if export_vars:
        journal = add_export_variables(journal, file_prefix, export_vars, export_zones)

    # Postprocess journal and check if there are any un-rendered variables
    journal = validate_journal(journal)

    return journal


def render_fluent_journal(template_path, template_map):
    """
    Read a template Fluent journal file and replace placeholders with values from the provided map.

    Read a template Fluent journal file and replace placeholders with values from the provided map.

    This function takes a template Fluent journal file where variables are specified in curly brackets `{}`.
    It then replaces these placeholders with the corresponding values provided in `template_map`.

    The function also checks that all keys from the `template_map` are used in the template. 
    If any keys are unused, the function raises a ValueError.

    Parameters
    ----------
    template_path : str
        Path to the template Fluent journal file containing placeholders in curly braces `{}`.
    template_map : dict
        Dictionary mapping the placeholders in the template to their replacement values. Values are expected to be numerical.

    Returns
    -------
    list of str
        List of strings where each string corresponds to a line from the rendered journal.

    Examples
    --------
    Suppose the content of 'journal_template.txt' is as follows:
    ```
    ; Define boundary conditions
    define/boundary-conditions/pressure-inlet inlet yes no {p0_in} no {p0_in} no yes no no yes {turbulence_intensity_in} {viscosity_ratio_in} 
    define/boundary-conditions/pressure-outlet outlet yes no {p_out} no yes no no yes {turbulence_intensity_in} {viscosity_ratio_in} no yes no no
    define/boundary-conditions/wall wall no no no no {wall_roughness_height} no {wall_roughness_constant}
    ```

    Using the function:

    >>> template_path = 'path_to_template.txt'
    >>> template_map = {
    ...     'p0_in': 2e6,
    ...     'p_out': 1e5,
    ...     'turbulence_intensity_in': 5,  # Percentage
    ...     'viscosity_ratio_in': 10, 
    ...     'wall_roughness_height': 1e-6,
    ...     'wall_roughness_constant': 0.5
    ... }
    >>> rendered_journal = render_fluent_journal(template_path, template_map)
    """
    
    # Read the template content
    with open(template_path, 'r') as file:
        template_content = file.readlines()

    # Convert all numerical values in template_map to strings for formatting
    string_template_map = {key: str(value) for key, value in template_map.items()}

    # Set of all keys in the template_map for tracking unused keys
    unused_keys = set(string_template_map.keys())

    rendered_journal = []
    for line in template_content:
        line = line.strip()  # Remove leading and trailing whitespaces
        if not line or line.startswith(';'):
            # If line is empty or a comment, add it as it is
            rendered_journal.append(line)
        else:
            # Check which keys are used in this line
            keys_used_in_line = [key for key in unused_keys if '{' + key + '}' in line]
            for key in keys_used_in_line:
                unused_keys.remove(key)
            # Replace content surrounded by curly brackets according to mapping
            rendered_journal.append(line.format(**string_template_map))

    # Check if there are any unused keys
    if unused_keys:
        raise ValueError(f"Some keys in the template_map were not used in the template: {', '.join(unused_keys)}.\nCheck for missing or mistyped variable names.")

    return rendered_journal


def add_barotropic_expressions(journal, expression_path, vars):
    """
    Modify the Fluent journal to add commands for defining barotropic model expressions.

    The function inserts the appropriate commands right after the "; Define Barotropic Model expressions" comment.

    Parameters:
    - journal (list of str): Existing journal entries.
    - expression_path (str): Path to the expression file.
    - vars (list of str): Variables to extract from the expression file.

    Returns:
    - list of str: Modified journal.

    Raises:
    - ValueError: If the "; Define Barotropic Model expressions" comment is not found.
    """

    # Remove the "barotropic_" prefix
    vars = [var.replace("barotropic_", "") for var in vars]

    # Read expressions from the provided file
    expressions = read_expression_file(expression_path, vars)
    
    # Create commands from expressions
    commands = []
    for var in vars:
        description = expressions[f"{var}_description"]
        expression = expressions[f"{var}_expression"]
        command = f'define/named-expressions edit "barotropic_{var}" description "{description}" definition "{expression}" q'
        commands.append(command)

    # Find the location to insert the commands in the journal
    try:
        idx = journal.index("; Define Barotropic Model expressions")
    except ValueError:
        raise ValueError('The "; Define Barotropic Model expressions" comment was not found in the journal.')

    # Insert the commands
    for i, command in enumerate(commands, start=idx+1):
        journal.insert(i, command)
    
    return journal




class read_expression_file:

    """
    Read content of barotropic model expression file, extracting desired variables' descriptions and expressions.

    The file should have:
    1. A header indicating the purpose of the file and its creation timestamp.
    2. Sections for each variable, with a description prefixed by ``case_`` followed by the mathematical expression.

    Parameters
    ----------
    filename : str
        Path to the barotropic model expression file.
    var_name : str
        List of variables to extract from the file.
    skiprows : int, optional
        Number of initial rows in the file to skip (default is 2 to bypass header).

    Returns
    -------
    str
        var_name.name = name of the variable. 
        var_name.expression = fluent expression of the variable.
    """

    def __init__(self, filename, var_name, skiprows=2):
        self.filename = filename
        self.var_name = var_name
        self.skiprow = skiprows

        # Check if the file exists
        if not os.path.isfile(self.filename):
            raise FileNotFoundError(f"The file '{self.filename}' does not exist. Stopping execution.")

    # Expression Description 
    def description(self):
        name = f"barotropic_model_{self.var_name}"
        return name



    def expression(self):
        with open(self.filename, "r") as file:
            lines = file.readlines()

        # Initialize a variable to store the expression as a string
        expression_str = ""
        start_collecting = False

        # Loop through each line, skipping header rows
        for line in lines[self.skiprow:]:
            # Check if the line starts with 'barotropic_model_{var_name}'
            if line.startswith(f"barotropic_model_{self.var_name}"):
                start_collecting = True  # Start collecting lines for the expression
                continue
            elif start_collecting:
                # Stop if we reach an empty line (end of expression)
                if line.strip() == "":
                    break
                # Append the line to the expression string
                expression_str += line.strip() + " "

        # return expression_str.strip()
        return expression_str



# def read_expression_file(filename, vars, skiprows=2):
#     """
#     Read content of barotropic model expression file, extracting desired variables' descriptions and expressions.

#     The file should have:
#     1. A header indicating the purpose of the file and its creation timestamp.
#     2. Sections for each variable, with a description prefixed by ``case_`` followed by the mathematical expression.

#     Parameters
#     ----------
#     filename : str
#         Path to the barotropic model expression file.
#     vars : list of str
#         List of variables to extract from the file.
#     skiprows : int, optional
#         Number of initial rows in the file to skip (default is 2 to bypass header).

#     Returns
#     -------
#     dict
#         Dictionary containing the description and expression for each desired variable. 
#         Each variable contributes two keys to the dictionary: '{var}_description' and '{var}_expression'.
#     """

#     # Open the file and read its lines
#     with open(filename, "r") as file:
#         lines = file.readlines()

#     # Initialize an empty dictionary to store the extracted expressions
#     expressions = {}
#     current_key = None

#     # Loop through each line, skipping header rows
#     for line in lines[skiprows:]:
#         line = line.strip()  # Remove leading and trailing whitespace

#         # Check if the line starts with 'case_', indicating it's a description line
#         if line.startswith("barotropic_"):
#             current_key = None  # Reset the current_key
            
#             # Determine if the description matches any of the desired variables
#             for var in vars:
#                 if var in line:
#                     current_key = var
#                     expressions[f"{var}_description"] = line
#                     expressions[f"{var}_expression"] = ""  # Initialize the expression as an empty string
#                     break

#         # If the line is part of an expression for a previously found variable, append it
#         elif current_key:
#             expressions[f"{current_key}_expression"] += line

#     return expressions


def add_solution_strategy(journal, strategy):
    """
    Insert the solution strategy into a Fluent journal after the designated comment.

    This function takes a Fluent journal list of commands and adds the solution strategy
    provided in the `strategy` dictionary right after the comment:
    `"; Run the simulation according to solution strategy"`.

    Parameters
    ----------
    journal : list of str
        List of strings, where each string corresponds to a command or comment in the initial Fluent journal content.
    strategy : dict
        A dictionary containing the solution strategy. It should have the following keys:

        - time_scale_factors : list of float
            List of time scale factors.
        - density_relaxation_factors : list of float
            List of density relaxation factors.
        - number_of_iterations : list of int
            List of number of iterations.

    Returns
    -------
    list of str
        List of strings with the modified Fluent journal commands.

    Raises
    ------
    ValueError
        If the "; Run the simulation according to solution strategy" comment is not found in the initial journal.

    Examples
    --------
    >>> journal_content = ["... [journal commands] ...", "; Run the simulation according to solution strategy", "... other commands ..."]
    >>> strategy = {
        "time_scale_factors": [0.25, 0.50],
        "density_relaxation_factors": [1.0, 1e-3],
        "number_of_iterations": [200, 300]
    }
    >>> new_journal = add_solution_strategy_to_journal(journal_content, strategy)
    """

    # Generate iteration commands as a list
    iteration_commands = []
    for timescale, relaxation, n_iter in zip(strategy["time_scale_factors"], 
                                             strategy["density_relaxation_factors"], 
                                             strategy["number_of_iterations"]):
        iteration_commands.append(f"solve/set/pseudo-time-method/global-time-step-settings yes 1 {timescale}")
        iteration_commands.append(f"solve/set/pseudo-time-method/relaxation-factors/density {relaxation}")
        iteration_commands.append(f"solve/iterate {n_iter}")

    # Find the index of the comment
    try:
        idx = journal.index("; Run the simulation according to solution strategy")
    except ValueError:
        raise ValueError('The "; Run the simulation according to solution strategy" comment was not found in the journal.')
    
    # Insert the generated iteration commands right after the comment
    for i, command in enumerate(iteration_commands, start=idx+1):
        journal.insert(i, command)

    return journal


def add_export_variables(journal, prefix, variable_mapping, zones):
    """
    Modify a Fluent journal to add commands for exporting specified variables to .xy files.

    The function takes an existing journal list and appends commands to export
    selected variables (as specified in variable_mapping) to .xy files. The filenames
    for the exported .xy files will be prefixed by the provided prefix and the variable name.
    The export commands are inserted right after the "; Save XY plot files" comment.

    Parameters
    ----------
    journal : list of str
        List of strings, where each string corresponds to a command or comment in the initial Fluent journal content.
    prefix : str
        Prefix for the exported .xy filenames, including the path.
    variable_mapping : dict
        Mapping from desired variable names (to be used in .xy filenames)
        to their corresponding names in Fluent.
    zones : list of str
        List of Fluent zone names or IDs to specify the zones of interest for exporting.

    Returns
    -------
    list of str
        List of strings with the modified Fluent journal commands.

    Raises
    ------
    ValueError
        If the "; Save XY plot files" comment is not found in the journal.

    Examples
    --------
    # Define mapping between export filenames and Fluent variable names:

    >>> EXPORT_VARS_MAPPING = {
            "pressure": "pressure",
            "density": "density",
            "viscosity": "viscosity-lam",
            "speed_sound": "sound-speed",
            ... # [Other mappings as needed]
        }

    # Define the journal content and prefix:

    >>> journal_content = ["... [journal commands] ...", "; Save XY plot files", "... other commands ..."]
    >>> prefix = "/path/to/simulation/case"

    # Add variables to export as .xy file to the journal:

    >>> new_journal = add_export_variables_to_journal(journal_content, prefix, EXPORT_VARS_MAPPING, ["axis", "wall"])
    """

    # Define commands to export .xy files
    export_entries = []
    for name, fluent_name in variable_mapping.items():
        filename = f"{prefix}_{name}.xy"
        export_entries.append(f"""plot/plot yes "{filename}" yes no no {fluent_name} yes 1 0 {' '.join(zones)} ()""")
        
    # Find the index of the comment
    try:
        idx = journal.index("; Save XY plot files")
    except ValueError:
        raise ValueError('The "; Save XY plot files" comment was not found in the journal.')
    
    # Insert the export commands right after the comment
    for i, command in enumerate(export_entries, start=idx+1):
        journal.insert(i, command)

    return journal


def validate_journal(journal):
    """
    Validate and combine the journal commands into a single string.

    This function takes a list of journal commands, combines them into a single
    string, and checks if there are any curly brackets left. If any curly brackets
    are found, it raises a ValueError indicating that not all variables were rendered.

    Parameters
    ----------
    journal : list of str
        List of strings where each string corresponds to a line from the processed journal.

    Returns
    -------
    str
        Single string containing all journal commands, separated by newlines.

    Raises
    ------
    ValueError
        If any curly brackets are found in the combined string, indicating that not all
        variables were successfully rendered.

    Examples
    --------
    >>> journal_lines = ["define bc pressure-inlet {p_in}", "solve"]
    >>> validate_and_combine_journal(journal_lines)
    ValueError: Not all variables were rendered in the journal. Check for missing or mistyped variable names.
    """

    # Filter out commented lines for validation
    non_commented_lines = [line for line in journal if not line.strip().startswith(';')]
    combined_journal = '\n'.join(non_commented_lines)

    # Find all occurrences of content within curly brackets
    unrendered_content = re.findall(r'\{([^}]+)\}', combined_journal) # chatgpt
    
    # Check if any curly brackets are left in the combined string
    if unrendered_content:
        error_message = "Not all variables were rendered in the journal. Unrendered variables found:\n"
        error_message += '\n'.join([f"  - {{{var}}}" for var in unrendered_content])
        raise ValueError(error_message)
    
    return '\n'.join(journal)



def parse_fluent_out(filename):
    """
    Parse an out-formatted datafile from Fluent and return a DataFrame.
    
    Parameters
    ----------
    filename : str
        The path to the Fluent output file.
        
    Returns
    -------
    pd.DataFrame
        A DataFrame containing the parsed data from the Fluent output file.
        
    Raises
    ------
    ValueError
        If the file format is not as expected.
    FileNotFoundError
        If the provided filename does not exist.
        
    Examples
    --------
    >>> df = parse_fluent_out("path_to_file.out")
    >>> print(df.head())
    """
    
    if not os.path.exists(filename):
        raise FileNotFoundError(f"The file {filename} was not found.")
    
    try:
        # Load file using Pandas
        df = pd.read_csv(filename, skiprows=2, delim_whitespace=True, header=0)

        # Process the column names
        cleaned_columns = [
            col.strip('"()').replace('(', '_').replace(')', '')
            for col in df.columns
        ]
        
        # Update the DataFrame columns and convert all values to floats
        df.columns = cleaned_columns
        df = df.apply(pd.to_numeric, errors='coerce')
        
    except Exception as e:
        raise ValueError(f"An error occurred while processing the file: {e}")
    
    return df


def parse_fluent_xy(filename):
    """
    Parse a Fluent .xy file and return its contents as a dictionary of DataFrames.

    Fluent's .xy files contain datasets, each representing XY plot data for some solution data 
    (like pressure or velocity) across a specific domain or region (like wall, axis, etc.). 
    This function reads such files and outputs a structured dictionary where the keys are dataset 
    labels (like "wall", "axis") and values are corresponding pandas DataFrames.

    Parameters
    ----------
    filename : str
        Path to the Fluent .xy file to be parsed.

    Returns
    -------
    data_out : dict
        A dictionary where:
        - The keys are dataset labels from the .xy file, e.g., "wall", "axis".
        - The values are pandas DataFrames containing the XY data. Each DataFrame has two columns 
        representing the X and Y axis data, labeled as per the .xy file's axis labels.

    Examples
    --------
    For a Fluent .xy file containing datasets "wall" and "axis", the output dictionary keys might be ['wall', 'axis'].
    """

    with open(filename, "r") as file:
        lines = file.readlines()

    # Extract axis labels
    result = parse('(labels "{x_label}" "{y_label}")\n', lines[1])
    x_label = result["x_label"]
    y_label = result["y_label"]

    data_out = {}
    zone_data = []
    zone_name = None

    # Iterate over the lines to extract datasets and their data points
    for line in lines[2:]:
        if line == "\n":  # skip empty lines
            continue
        
        # Detect start of a new dataset
        if line.startswith("((xy/key/label"):
            # If previous data exists, store it in the dictionary before starting a new dataset
            if zone_name and zone_data:
                data_out[zone_name] = pd.DataFrame(zone_data, columns=[x_label, y_label])
                zone_data = []

            zone_name = parse('((xy/key/label "{zone}")\n', line)["zone"]
        
        # Detect end of current dataset
        elif line.startswith(")"):
            if zone_name and zone_data:
                data_out[zone_name] = pd.DataFrame(zone_data, columns=[x_label, y_label])
                zone_data = []
        
        # Extract XY data points
        else:
            x, y = map(float, line.split('\t'))  # Conver all values to numeric
            zone_data.append([x, y])

    return data_out


def plot_residuals(transcript_file, fig=None, ax=None, savefig=False, fullpath=None):
    """
    Plot the residual history from the given transcript file.

    Parameters
    ----------
    transcript_file : str
        Path to the transcript file containing residuals.
    fig : matplotlib.figure.Figure, optional
        Existing figure to plot on, if provided. Defaults to None.
    ax : matplotlib.axes.Axes, optional
        Existing axis to plot on, if provided. Defaults to None.
    savefig : bool, optional
        Flag to save the figure to a file. Defaults to False.
    fullpath : str, optional
        If savefig is True, path to save the figure to (excluding file extension). Defaults to None.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure containing the plotted residuals.
    ax : matplotlib.axes.Axes
        Axes containing the plotted residuals.

    Notes
    -----
    The function reads the residuals from the given transcript file using `read_residual_file` 
    and then plots them on a semi-log scale. If `fig` and `ax` are not provided, 
    it creates a new figure and axis. The residuals are mapped to readable labels 
    using a predefined name_map.
    """

    # Read residuals from the transcript file
    df = read_residual_file(transcript_file, res_number=5)

    # If no figure or axis provided, create new ones
    if not fig and not ax:
        fig, ax = plt.subplots(figsize=(6.4, 4.8))

    # Set plot title, labels and scales
    ax.set_title(f'Convergence history from file {os.path.basename(transcript_file)}')
    ax.set_xlabel("Iteration number")
    ax.set_ylabel("Normalized RMS residual")
    ax.set_yscale("log")

    # Define a mapping for the residual names to readable labels
    # TODO: extend to 3D and other turbulence models
    name_map = {
        'continuity': '$p$ - continuity',
        'x-velocity': '$v_x$ - momentum',
        'y-velocity': '$v_y$ - momentum',
        'k': '$k$ - turbulence',
        'omega': '$\omega$ - turbulence'
    }

    # Plot residuals using the mapped labels
    for column in df.columns[1:]:
        ax.plot(
            df["iter"], df[column],
            linewidth=1.0, linestyle="-", marker="none", 
            label=name_map[column]
        )

    # Add legend to the plot
    ax.legend(loc="upper right", fontsize=11)

    # Adjust plot layout
    fig.tight_layout(pad=1)

    # If save flag is set and a path is provided, save the figure in multiple formats
    if savefig and fullpath:
        for ext in [".png", ".svg", ".pdf"]:
            fig.savefig(fullpath + ext, bbox_inches="tight")

    return fig, ax


def read_residual_file(filename, res_number=5, line_start=0, header=None):
    """
    Reads a transcript file to extract the residual history and returns it as a Pandas DataFrame.

    Parameters
    ----------
    filename : str
        Path to the transcript file.
    res_number : int, optional
        Number of residuals to be extracted. Defaults to 5.
    line_start : int, optional
        Line number to start reading from. Useful for files with non-relevant introductory lines. Defaults to 0.
    header : list of str, optional
        List of column names for the DataFrame. If not provided, the function tries to extract it from the file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the extracted residuals. Column names are either provided or extracted from the file.

    Notes
    -----
    It's assumed that the residual headers in the transcript file start with 'iter' and 'continuity'.
    """
    
    # Read data from the transcript file
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Initialize lists to store values
    residuals = []
    
    # Use the provided header if available, otherwise extract from the file
    if not header:
        header = []
    
    # Loop through lines to extract headers and residuals
    for line in lines[line_start:]:
        
        # If header hasn't been found, look for it based on the expected starting words
        words = line.split()
        if not header and len(words) >= 2 and words[0] == 'iter' and words[1] == 'continuity':
            header = line.split()
            header = [item.split("/") for item in header]
            header = [item for sublist in header for item in sublist]
            header = header[0:res_number+1]

        # Extract lines that contain numeric residual values
        line_elements = line.split()[0:res_number+1]
        numeric_line = all([is_float(element) for element in line_elements])
        
        # Append lines that match the expected structure and content
        if header and len(header) == len(line_elements) and numeric_line:
            residuals.append(line_elements)

    # Convert the list of residuals into a DataFrame
    df = pd.DataFrame(residuals, columns=header)

    # Convert values to float for computation purposes
    df = df.astype(float)

    return df


def print_transcript_real_time(transcript_file):
    """Print the content of transcript_file.trn in real time"""

    # Wait until the transcript file is created before executing the script
    wait_for_file(transcript_file)

    # Run a subprocess to run the "stream_tramscript_file.py" script in the background
    args = ["python", stream_transcript.__file__, f"{transcript_file}"]
    process = subprocess.Popen(args)
    
    return process
    

def print_transcript_real_time(transcript_file):
    """
    Print the content of transcript_file.trn in real time.

    This function launches the `stream_transcript.py` script as a subprocess 
    to print the content of the given transcript file in real time.

    Parameters
    ----------
    transcript_file : str
        Path to the transcript file to be printed.

    Returns
    -------
    subprocess.Popen
        The subprocess object representing the background process.
    """
    wait_for_file(transcript_file)
    args = [sys.executable, stream_transcript.__file__,  transcript_file]
    return subprocess.Popen(args)

def plot_residuals_real_time(transcript_file):
    """
    Plot the residual history contained in the transcript_file.trn in real time.

    This function launches the `stream_residuals.py` script as a subprocess 
    to plot the residual history of the given transcript file in real time.

    Parameters
    ----------
    transcript_file : str
        Path to the transcript file whose residuals are to be plotted.

    Returns
    -------
    subprocess.Popen
        The subprocess object representing the background process.
    """
    wait_for_file(transcript_file)
    args = [sys.executable, stream_residuals.__file__,  transcript_file]
    return subprocess.Popen(args)


def run_fluent_simulation(journal, n_proc=1, timeout=3600, plot_residuals=False, 
                          print_transcript=False, transcript_file=None):
    """
    Orchestrates the execution of a Fluent simulation using the commands from the provided journal file.
    
    The function manages the entire simulation process from initiating Fluent, sending it the right 
    commands, visualizing output, to cleaning up post-simulation resources. It offers functionalities 
    to both visualize the Fluent console output in Python and plot the residual history in real time.

    Parameters
    ----------
    journal : str
        The sequence of Fluent commands to be executed.
    n_proc : int, optional
        Number of processor cores Fluent should use. Default is 1.
    timeout : int, optional
        Maximum allowable time for the Fluent process to run, in seconds. Default is 3600 (1 hour).
    plot_residuals : bool, optional
        If True, plots residuals in real time. Requires `transcript_file` to be provided. Default is False.
    print_transcript : bool, optional
        If True, the Fluent console output will be displayed in Python. Requires `transcript_file` to be provided. Default is False.
    transcript_file : str, optional
        Path to the transcript file which contains the residuals and other output data from Fluent. 
        Mandatory if either `plot_residuals` or `print_transcript` are set to True.

    Raises
    ------
    TranscriptFileMissingError
        If either `plot_residuals` or `print_transcript` is True and `transcript_file` is not provided.

    Returns
    -------
    None

    Notes
    -----
    After the simulation, the function automatically invokes the Fluent cleanup script to close open processes.
    """
    

    # Check for FLUENT_EXEC environment variable
    fluent_exec_path = os.getenv("FLUENT_EXEC")
    if fluent_exec_path is None:
        raise EnvironmentError(
            "Environment variable FLUENT_EXEC is not set. Please add the following lines to your .bashrc file (or your shell's equivalent), "
            "adjusting the paths to match the location where Fluent is installed on your system:\n\n"
            "# Path to the Fluent v23.1 bin directory\n"
            "export FLUENT_PATH=\"/path/to/your/ansys/inc/v231/fluent/ntbin/win64\"\n"
            "export FLUENT_EXEC=\"$FLUENT_PATH/fluent.exe\"\n"
            "export PATH=$PATH:$FLUENT_PATH\n\n"
            "Then restart your shell or run 'source ~/.bashrc' (without quotes) for the changes to take effect."
        )

    # Check if the Fluent executable exists at the specified path
    if not shutil.which(fluent_exec_path):
        raise FileNotFoundError(
            f"Fluent executable not found at the path specified by the FLUENT_EXEC environment variable: {fluent_exec_path}. "
            "Please ensure that the FLUENT_EXEC environment variable is set correctly in your .bashrc file (or your shell's equivalent), "
            "and points to the Fluent executable."
        )

    plot_process = None
    print_process = None
    exit_flag = None
    cleanup_done = False

    # Run simulations
    try:

        # Create Fluent child process
        # 2dpp -> 2D solver with double precision arithmetic
        # -g -> Do not start up the graphical user interface
        # -tN -> Run the simulation using N processor cores
        process = pexpect.popen_spawn.PopenSpawn(
            f'{os.getenv("FLUENT_EXEC")} 2ddp -g -t{n_proc}'
        )

        # Send sequence of commands to Fluent CLI
        process.send(journal)

        # Plot residuals from transcript file in real time
        if plot_residuals:
            if transcript_file is None:
                raise TranscriptFileMissingError()
            plot_process = plot_residuals_real_time(transcript_file)  # Function waits until the transcript file exists

        # Print the transcript file in real time
        if print_transcript:
            if transcript_file is None:
                raise TranscriptFileMissingError()
            print_process = print_transcript_real_time(transcript_file)  # Function waits until the transcript file exists

        # Wait for the child process to finish
        exit_flag = process.expect([pexpect.EOF, pexpect.TIMEOUT], timeout=timeout)
        # print(process.before.decode('utf-8'))

        # Check execution status
        if exit_flag == 0:
            # Fluent execution reached the End Of File (EOF)
            print("Simulation completed\n")

        elif exit_flag == 1:
            # Time out reached before Fluent execution was complete
            print("Timeout before execution was complete")

    except KeyboardInterrupt:
        print("Fluent execution terminated by user.")
        if not cleanup_done:
            run_cleanup_script()
            cleanup_done = True
        print("Continuing script execution...")
        return # Exit the function

    except TranscriptFileMissingError as e:
        if not cleanup_done:
            run_cleanup_script()
            cleanup_done = True
        raise e

    finally:
        # Cleanup operations
        if plot_process:
            plot_process.terminate()

        if print_process:
            print_process.terminate()

        if exit_flag != 0 and not cleanup_done:
            run_cleanup_script()

 
class TranscriptFileMissingError(ValueError):
    """Raised when plotting residuals is requested but no transcript file is provided."""
    
    def __init__(self, message=None):
        if message is None:
            message = ("To plot residuals in real-time, a transcript file is required. "
                       "Please provide the 'transcript_file' parameter when setting 'plot_residuals' to True.")
        super().__init__(message)


def run_cleanup_script(max_retries=5, wait_interval=5):
    """
    Execute the cleanup script found in the current working directory.

    This function repeatedly attempts to find and execute files with a prefix 
    of "cleanup-fluent" in the current directory. If the script isn't found 
    immediately, the function will wait for a short interval and try again, 
    up to a specified maximum number of retries.

    Parameters:
        max_retries (int): Maximum number of times the function will attempt 
                           to find and execute the cleanup script.
        wait_interval (int): Time in seconds to wait between consecutive 
                             retries.

    Note:
        The function is used to ensure that necessary cleanup operations are 
        performed after the simulation, especially when the creation of the 
        cleanup script might be delayed.
    """

    print("Fluent did not terminate properly. Using the Fluent cleanup script...")
    CLEANUP_PREFIX = "cleanup-fluent"
    retries = 0

    while retries < max_retries:
        for filename in os.listdir(os.getcwd()):
            if filename.startswith(CLEANUP_PREFIX):
                os.system(filename)
                return
        retries += 1
        time.sleep(wait_interval)  # Wait before the next retry

    print("Failed to find cleanup script after", max_retries, "retries.")










