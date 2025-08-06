import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import inspect

def plotsettings():
    font = 22
    plt.rcParams['font.size'] = font  # Default font size (affects labels, titles, legend)
    plt.rcParams['lines.linewidth'] = 2  # Default line width
    plt.rcParams['lines.markersize'] = 10  # Default marker size
    plt.rcParams['legend.fontsize'] = font  # Default legend font size
    plt.rcParams['xtick.labelsize'] = font  # X-axis tick label size
    plt.rcParams['ytick.labelsize'] = font  # Y-axis tick label size
    mpl.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['figure.figsize'] = (10,8)

def gotoFilePath():
    #current_file_path = os.path.abspath(__file__)
    #workdir = os.path.dirname(current_file_path)
    #os.chdir(workdir)

    # Get the file of the caller
    caller_frame = inspect.stack()[1]
    caller_path = os.path.abspath(caller_frame.filename)
    caller_dir = os.path.dirname(caller_path)
    os.chdir(caller_dir)    
