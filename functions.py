import os
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from getpass import getuser
from numba import njit

@njit
def time_series(x0: np.float64, y0: np.float64, a: np.float64, b: np.float64, c: np.float64, m: np.float64, N: np.int32) -> np.ndarray:
    
    x = x0
    y = y0
    u = np.zeros((N, 2))


    for i in range(N):
        y = y - b * np.sin(2 * np.pi * x) - c * np.sin(2 * np.pi * m * x)
        x = (x + a * (1 - y ** 2)) % 1.0
        u[i, 0] = x
        u[i, 1] = y

    return u

def time_series_ICs(x0: np.ndarray, y0:np.ndarray, a: np.float64, b: np.float64, c: np.float64, m: np.float64, N: np.int32):

    num_ic = x0.shape[0]
    u = np.zeros((num_ic * N, 2))
    for i in range(num_ic):
        u[i * N:(i + 1) * N] = time_series(x0[i], y0[i], a, b, c, m, N)
    
    return u


def plot_params(fontsize=20, legend_fontsize=14, axes_linewidth=1.3):
    """
    Update the parameters of the plot.

    Returns
    -------
    cmap : string
        The color map used in the colored plots.
    """
    tick_labelsize = fontsize - 3

    plt.clf()
    plt.rc('font', size=fontsize)
    plt.rc('xtick', labelsize=tick_labelsize)
    plt.rc('ytick', labelsize=tick_labelsize)
    plt.rc('legend', fontsize=legend_fontsize)
    font = {'family' : 'STIXGeneral'}
    plt.rc('font', **font)
    plt.rcParams["mathtext.fontset"] = "stix"
    mpl.rcParams['axes.linewidth'] = axes_linewidth #set the value globally

def set_ticks_in(ax, pad_x = False, pad_y = False):

    try:
        if len(ax.shape) == 1:
            for i in range(ax.shape[0]):
                ax[i].tick_params(axis="both", direction='in', which="major")
                ax[i].tick_params(axis="both", direction='in', which="minor")
                if pad_x:
                    ax[i].tick_params(axis="x", which="major", pad=pad_x)
                if pad_y:
                    ax[i].tick_params(axis="y", which="major", pad=pad_y)
        else:
            for i in range(ax.shape[0]):
                for j in range(ax.shape[1]):
                    ax[i, j].tick_params(axis="both", direction='in')
                    if pad_x:
                        ax[i, j].tick_params(axis="x", which="major", pad=pad_x)
                    if pad_y:
                        ax[i, j].tick_params(axis="y", which="major", pad=pad_y)
    except  AttributeError:
        ax.tick_params(axis="both", direction='in', which="major")
        ax.tick_params(axis="both", direction='in', which="minor")
        if pad_x:
            ax.tick_params(axis="x", which="major", pad=pad_x)
        if pad_y:
            ax.tick_params(axis="y", which="major", pad=pad_y)

def truncate_colormap(cmap="nipy_spectral", minval=0.0, maxval=0.95, n=256):
    """Truncate a colormap to only include a subset of its range."""
    original_cmap = plt.get_cmap(cmap)
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        f'trunc({original_cmap.name},{minval:.2f},{maxval:.2f})',
        original_cmap(np.linspace(minval, maxval, n))
    )
    return new_cmap

def get_path(current_dir, subdirs=None):
    user = getuser()
    
    # Determine the base path based on the user
    if user in {"jdanilo", "jdsjunior"}:
        base_path = os.path.join("/Users", user, "Matheus", "Pesquisa", current_dir)
    else:
        base_path = os.path.join("/Users", user, "Pesquisa", current_dir)
    
    # Append subdirs if provided
    if subdirs:
        path = os.path.join(base_path, subdirs)
    else:
        path = base_path
    
    # Create the directory if it doesn't exist
    if not os.path.isdir(path):
        print(f"\t{path} does not exist! Trying to create it...")
        os.makedirs(path, exist_ok=True)
        if not os.path.isdir(path):
            print(f"\t{path} could not be created. Exiting...")
            sys.exit(1)
    
    return path

def extract_grid(datafile):
    df = pd.read_csv(datafile, header=None, sep=" ")
    x = np.array(df[0])
    y = np.array(df[1])
    z = np.array(df[2])
    grid = int(np.sqrt(len(x)))
    x = x.reshape((grid, grid))
    y = y.reshape((grid, grid))
    z = z.reshape((grid, grid))

    return x, y, z

def generate_sbatch_script(input_files, output_files, script_name, args, time="24:00:00", email="rolim.sales@unesp.br", start=False, end=False, clusters=1):
    """
    Generate an sbatch script as a string.
    
    Parameters:
        email (str): Email address for job notifications.
        input_files (list): List of input files.
        output_files (list): List of output file patterns.
        script_name (str): Name of the Python script to run.
        args (list): List of arguments to pass to the Python script.
        time (str): Wall time in SLURM format (default: "4-00").
        
    Returns:
        str: The sbatch script as a string.
    """

    if not isinstance(clusters, int) or clusters < 1:
        raise ValueError("clusters must be a positive integer.")
   
    input_files_str = " ".join(input_files)
    output_files_str = " ".join(output_files)
    args_str = " ".join(f"${{args[{i}]}}" for i in range(len(args)))
    
    if script_name.endswith(".py"):
        comm = "python "
        export_statements="""
export PATH=/home/rolim/anaconda3/bin:$PATH
export LD_LIBRARY_PATH=/home/rolim/anaconda3/lib:$LD_LIBRARY_PATH
export INCLUDE=/home/rolim/anaconda3/include:$INCLUDE
"""
    elif script_name.endswith(".x"):
        comm = "./"
        export_statements=""
    else:
        raise ValueError("Invalid script extension. Only .py and .x are supported.")
    
    if clusters == 1:
        time_statement = f"#SBATCH -t {time}"
    else:
        time_statement = f"#SBATCH -t {time} -c {clusters}"
    
    if email != False:
        if start != False and end != False:
            email_statement = f"#SBATCH --mail-user={email} --mail-type=ALL"
        elif start != False:
            email_statement = f"#SBATCH --mail-user={email} --mail-type=BEGIN"
        elif end != False:
            email_statement = f"#SBATCH --mail-user={email} --mail-type=END,FAIL"
        else:
            email_statement = ""
    else:
        email_statement = ""

    script = f"""#!/bin/bash

{time_statement}
{email_statement}

export INPUT="{input_files_str}"
export OUTPUT="{output_files_str}"

{export_statements}

args=("$@")

job-nanny {comm}{script_name} {args_str}
"""
    return script

def format_number(x, decimal_places=5):
    formatted = f"{x:.{decimal_places}f}"
    if formatted.startswith('0'):
        formatted = formatted[1:]
    elif formatted.startswith("-0"):
        formatted = formatted[0] + formatted[2:]
    return formatted

def should_process(i):

    process = False

    if i <= 1e4:
        if (i - 1) % 2 == 0:
            process = True
    elif 1e4 < i <= 1e5:
        if i % 500 == 0:
            process = True
    elif 1e5 < i <= 1e6:
        if i % 5000 == 0:
            process = True
    elif 1e6 < i <= 1e7:
        if i % 50000 == 0:
            process = True
    elif 1e7 < i <= 1e8:
        if i % 5e5 == 0:
            process = True
    elif 1e8 < i <= 1e9:
        if i % 5e6 == 0:
            process = True
    elif 1e9 < i <= 1e10:
        if i % 5e7 == 0:
            process = True
    elif 1e10 < i <= 1e11:
        if i % 5e8 == 0:
            process = True
    elif 1e11 < i <= 1e12:
        if i % 5e9 == 0:
            process = True

    return process

def generate_values_around_integer(n):
    """
    Generates a list of values around the given integer n.
    
    Parameters:
        n (int): The integer around which to generate values.
    
    Returns:
        list: A list of values around the integer n.
    """
    # Generate values below n
    below_n = np.arange(n - 0.05, n, 0.01).tolist() + np.arange(n - 0.005, n, 0.001).tolist()
    
    # Generate values above n
    above_n = [n] + np.arange(n + 0.001, n + 0.006, 0.001).tolist() + np.arange(n + 0.01, n + 0.06, 0.01).tolist()
    
    # Combine the lists
    values = below_n + above_n
    
    # Sort the list to ensure the correct order
    values.sort()
    
    return values