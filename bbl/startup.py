import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime
import random
import numpy as np

font = 20
font = 18
plt.rcParams['font.size'] = font
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 10
plt.rcParams['legend.fontsize'] = font
plt.rcParams['xtick.labelsize'] = font
plt.rcParams['ytick.labelsize'] = font
mpl.rcParams['savefig.bbox'] = 'tight'

# plt.rcParams['figure.figsize'] = (10,8)
plt.rcParams['figure.figsize'] = (8,6)

cla = lambda: plt.close('all')


# # auto switch into the directory of the script
# # ============================================== 
# import os
# os.chdir(os.path.dirname(os.path.abspath(__file__)))


def clearall():
    """Fully reset the IPython interactive namespace (like MATLAB 'clear all')"""
    from IPython import get_ipython
    ipy = get_ipython()
    if ipy:
        ipy.run_line_magic('reset', '-f')
        
        
import subprocess
def winpath(path):
    if '\x00' in path:
        raise ValueError("Path contains null byte (\\x00), which is not allowed.")

    try:
        result = subprocess.run(
            ['wslpath', '-a', path],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print("wslpath failed:", e.stderr)
        return None
# win_path = r"C:\Users\Username\Documents\file.txt"

# def setplot(figsize=(12,10)):
def setplot(figsize=(8,6), font=18):
    # font = 20
    # font = 18
    plt.rcParams['font.size'] = font
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['lines.markersize'] = 10
    plt.rcParams['legend.fontsize'] = font
    plt.rcParams['xtick.labelsize'] = font
    plt.rcParams['ytick.labelsize'] = font    
    mpl.rcParams['savefig.bbox'] = 'tight'

    # plt.rcParams['figure.figsize'] = (10,8)
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['axes.grid'] = True
    
def savefig(filename=None):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    jj = random.randint(0,10)
    
    if filename == None:
        plt.savefig(f"{timestamp}_{jj}.png", dpi=300, bbox_inches='tight')    
    else:
        plt.savefig(f"{filename}_{timestamp}_{jj}.png", dpi=300, bbox_inches='tight')     
    
from matplotlib.ticker import LogLocator, LogFormatter    
def set_logy_ticks(ax):
    #add minor ticks on y axis 
    # ax.yaxis.set_major_locator(LogLocator(base=10.0))
    # ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2,10)*0.1, numticks=10))
    # ax.yaxis.set_minor_formatter(LogFormatter(base=10.0, labelOnlyBase=False))    
    
    # ax = plt.gca()
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2,10), numticks=15))
    # ax.tick_params(which='major', length=8)
    # ax.tick_params(which='minor', length=4)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    
    
   