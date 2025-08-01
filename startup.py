import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime
import random

font = 20
font = 15
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

def setplot():
    font = 20
    font = 15
    plt.rcParams['font.size'] = font
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['lines.markersize'] = 10
    plt.rcParams['legend.fontsize'] = font
    plt.rcParams['xtick.labelsize'] = font
    plt.rcParams['ytick.labelsize'] = font    
    mpl.rcParams['savefig.bbox'] = 'tight'

    # plt.rcParams['figure.figsize'] = (10,8)
    plt.rcParams['figure.figsize'] = (8,6)
    
    
def savefig():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    jj = random.randint(0,10)

    plt.savefig(f"{timestamp}_{jj}.png", dpi=300, bbox_inches='tight')    
    
    
   