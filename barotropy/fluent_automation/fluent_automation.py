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


    
def read_expression_file(filename):
    
    '''
    Read all the expressions contained in the output file of the barotropic model and store them in a dictionary with the name of the variable as "barotropic_{variable_name}"
    The expression file is needed as input
    '''

    # Read the file and store all in lines
    with open(filename, "r") as file:
        lines = file.readlines()

    # Initialize the dictionary
    expressions = {}
    start_collecting = False

    for line in lines:

        if line.startswith("barotropic_"):
            name = line.strip()
            start_collecting = True
            expressions[name] = []
            expression_str = ""
            continue

        elif start_collecting:
                if line.strip() == "":
                    start_collecting = False
                    expressions[name] = expression_str.strip()
                    continue
                # Append the line to the expression string
                expression_str += line.strip() + " "

    return expressions


    


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
        'continuity': r'$p$ - continuity',
        'x-velocity': r'$v_x$ - momentum',
        'y-velocity': r'$v_y$ - momentum',
        'k': r'$k$ - turbulence',
        'omega': r'$\omega$ - turbulence'
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


# def print_transcript_real_time(transcript_file):
#     """Print the content of transcript_file.trn in real time"""

#     # Wait until the transcript file is created before executing the script
#     wait_for_file(transcript_file)

#     # Run a subprocess to run the "stream_tramscript_file.py" script in the background
#     args = ["python", stream_transcript.__file__, f"{transcript_file}"]
#     process = subprocess.Popen(args)
    
#     return process
    

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










