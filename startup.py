import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

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

def clearall():
    """Fully reset the IPython interactive namespace (like MATLAB 'clear all')"""
    from IPython import get_ipython
    ipy = get_ipython()
    if ipy:
        ipy.run_line_magic('reset', '-f')
