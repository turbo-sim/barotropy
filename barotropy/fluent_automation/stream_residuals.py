import os
import sys
import time
import pandas as pd

import matplotlib.pyplot as plt


RESIDUAL_LABEL_MAPPING = {
    'continuity': r'$p$ - continuity',
    'x-velocity': r'$v_x$ - momentum',
    'y-velocity': r'$v_y$ - momentum',
    'z-velocity': r'$v_z$ - momentum',
    'k': r'$k$ - turbulence',
    'omega': r'$\omega$ - turbulence'
}


def _plot_residuals_real_time(filename, frequency=0.5):
    """Plot the residual history from transcript file in real time"""
    print("Plotting data")

    # Import within function to prevent circular import
    from barotropy import read_residual_file, set_plot_options

    # Set nice plot settings
    set_plot_options(fontsize=14, grid=True, margin=0.00)

    # Read the transcript file until there is residual information
    df = pd.DataFrame()
    while df.empty:
        df,_ = read_residual_file(filename)
        time.sleep(1)  # Adjust the delay time as needed
    
    # Create the figure and initialize lines for each residual
    fig, ax = plt.subplots()
    lines = {}
    for res_name in df.columns[1:]:
        line, = ax.plot([], [], linewidth=1.0, label=RESIDUAL_LABEL_MAPPING.get(res_name, res_name))
        lines[res_name] = line

    # Define title and axes labesl
    ax.set_title('Convergence history from file: ' + os.path.basename(filename))
    ax.set_xlabel("Iteration number")
    ax.set_ylabel("Normalized RMS residual")
    ax.set_xscale("linear")
    ax.set_yscale("log")
    
    # Create legend
    ax.legend(loc="upper right", fontsize=11)

    # Adjust PAD
    fig.tight_layout(pad=1, w_pad=None, h_pad=None)

    # Prepare interactive figure
    plt.ion()  # Turn on interactive mode
    plt.show()  # Display the initial empty plot

    try:
        while plt.get_fignums(): # Continue loop while figures are open

            # Read residual history from file
            df,_ = read_residual_file(filename)

            # Update the plot with the latest data
            for res_name, line in lines.items():
                line.set_xdata(df["iter"].values)
                line.set_ydata(df[res_name].values)

            # Update the axes limits
            ax.relim()
            ax.autoscale_view()

            # Draw the figure
            plt.draw()
            plt.pause(frequency)  # Add a short pause to allow the plot to update


    except Exception as e:
        print(f"Error: {e}")

    except KeyboardInterrupt:
        pass

    finally:
        plt.close()  


if __name__ == "__main__":

    # Get the transcript file as command-line argument
    if len(sys.argv) != 2:
        print("Usage: python main_script.py <file_path>")
        sys.exit(1)
    
    # Plot the residual history as the file gets updated
    _plot_residuals_real_time(filename=sys.argv[1])