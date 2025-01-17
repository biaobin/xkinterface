import sys
#sys.path.append('/afs/ifh.de/group/pitz/data/lixiangk/work/sync/python')

#from astra_modules import *
from interface import *
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#plt.switch_backend('agg')

'''
General plot functions for Astra and SCO.
Options:
  ev: evolutions of beam properties stored in xxx.Xemit.run
  ps: phase space of a given beam distribution
  sco: evolutions of beam properties stored in BeamDynamics.dat
Examples:
`
python $Py/demo1.py sco suffix='1__.txt' fig_ext='_1.eps'
`
'''

def plot_config():
    from cycler import cycler
    from matplotlib.ticker import AutoMinorLocator

    fsize = 14 # a quarter of the paper width: 20 pt; half of the paper width: 12
    font = {'size' : fsize, 'family' : 'serif'}
    color_cycle = ['r', 'b', 'g', 'c', 'm', 'y', 'k']
    linestyle_cycle = ['-', '--', '-.', ':', (0, (5, 2, 5, 2)), (0, (10, 2, 5, 2, 2, 2)), (0, (12, 2, 2, 2))]
    marker_cycle = ['o', 'd', 'v', '^', '<', '>', '*']
    
    mpl.rc('font', **font)
    mpl.rc('xtick', labelsize = 10, direction = 'in', top   = True)
    mpl.rc('ytick', labelsize = 10, direction = 'in', right = True)
    mpl.rc('xtick.major', size = 5, width = 1)
    mpl.rc('ytick.major', size = 5, width = 1)
    mpl.rc('xtick.minor', size = 3, width = 0.7, visible = True)
    mpl.rc('ytick.minor', size = 3, width = 0.7, visible = True)
    
    mpl.rc('lines', linewidth=2, markersize=6, color='r')
    # mpl.rc('lines', linestyle = 'solid')
    mpl.rc('axes', labelpad = 0, prop_cycle=(cycler('color', color_cycle) + cycler('linestyle', linestyle_cycle) + cycler('marker', marker_cycle)))
    mpl.rc('legend', fontsize = 12, labelspacing = 0.05, handletextpad=0.4, frameon=False, handlelength=2.1)
    
    mpl.rc('figure', dpi = 100, figsize = (4, 4))
    mpl.rc('figure.subplot', bottom = 0.15, top = 0.9, left = 0.15, right = 0.9)
    
    mpl.rc('image', cmap = 'jet')
    
    return
plot_config()


argv = sys.argv
print(argv)
kv = dict([arg.split('=', 1) for arg in argv[2:]])
print(kv)

if 'workdir' in list(kv.keys()):
    workdir = kv['workdir']
else:
    workdir = r'./'
    
if 'fig_ext' in list(kv.keys()):
    fig_ext = kv['fig_ext']
else:
    fig_ext = '.eps'

Lextent = True
if 'xmax' in list(kv.keys()):
    xmax = np.float(kv['xmax'])
    if 'xmin' in list(kv.keys()):
        xmin = np.float(kv['xmin'])
    else:
        xmin = -xmax
else:
    Lextent = False
if 'ymax' in list(kv.keys()):
    ymax = np.float(kv['ymax'])
    if 'ymin' in list(kv.keys()):
        ymin = np.float(kv['ymin'])
    else:
        ymin = -ymax
else:
    Lextent = False

if not Lextent:
    extent = None
else:
    extent = [xmin, xmax, ymin, ymax]
print(extent)

os.chdir(workdir)

if argv[1] in ['evolution', 'ev']:

    if 'prefix' in list(kv.keys()):
        prefix = kv['prefix']
    else:
        prefix = 'ast'

    if 'suffix' in list(kv.keys()):
        suffix = kv['suffix']
    else:
        suffix = '001'
    fig_ext = '@%s.%s' % (suffix, fig_ext)
    
    #plot_avg_xy(prefix = prefix, suffix = suffix, extent = None, figsize=(8, 4))
    plot_rms_xy(prefix = prefix, suffix = suffix, fig_ext = fig_ext, extent = extent, figsize=(8, 4))
    plot_rms_xyz(prefix = prefix, suffix = suffix, fig_ext = fig_ext, extent = extent, figsize=(8, 4))
    plot_emi_xy(prefix = prefix, suffix = suffix, fig_ext = fig_ext, extent = None, figsize=(8, 4))
    plot_kin(   prefix = prefix, suffix = suffix, fig_ext = fig_ext, extent = None, figsize=(8, 4))
    #sys.exit()

if argv[1] in ['phase space', 'ps']:

    if 'fname' in list(kv.keys()):
        fname = kv['fname']
    else:
        fname = 'ast.0528.001'
        
    if 'vmin' in list(kv.keys()):
        vmin = kv['vmin']
    else:
        vmin = 0.005
        
    fig_ext = '@%s%s' % (fname, fig_ext)
    
    diag = BeamDiagnostics(fname = fname, plot = True, fig_ext = fig_ext, energy = False)
    diag.demo(fname)
    
    xmax = ymax = 4
    plot2d_xy( dist = diag.dist, fig_ext = fig_ext, vmin = vmin, bins = (200, 200), extent = extent)
    plot2d_xpx(dist = diag.dist, fig_ext = fig_ext, vmin = vmin, bins = (200, 200))
    plot2d_ypy(dist = diag.dist, fig_ext = fig_ext, vmin = vmin, bins = (200, 200))
    plot2d_zpz(fname, fig_ext = fig_ext, vmin = 0.01, bins = (200, 200))
    #sys.exit()

if argv[1] in ['sco']:

    if 'fname' in list(kv.keys()):
        fname = kv['fname']
    else:
        fname = 'BeamDynamics'
        
    if 'suffix' in list(kv.keys()):
        suffix = kv['suffix']
    else:
        suffix = '.dat'
    
    if 'z0' in list(kv.keys()):
        z0 = kv['z0']
        z0 = float(z0)
    else:
        z0 = 0
    
    fname = fname+suffix

    plot_sco(fname, extent = extent, figsize=(8, 4), figname = 'sco-rms_xy-z', z0 = z0, fig_ext = fig_ext)
    plot_sco_xyz(fname, extent = extent, figsize=(8, 4), figname = 'sco-rms_xyz-z', z0 = z0, fig_ext = fig_ext)
    plot_sco(fname, y = ['nemit_x', 'nemit_y'], ylabel = u_emi_xy, figname = 'sco-emi_xy-z', z0 = z0, fig_ext = fig_ext,\
             extent = None, figsize=(8, 4))
    #sys.exit()

if argv[1] in ['slice', 'sl']:

    if 'fname' in list(kv.keys()):
        fname = kv['fname']
    else:
        fname = 'ast.0528.001'
    fig_ext = '@%s.eps' % fname
    
    d1 = astra2slice(fname, 'slice@'+fname, nslice = 20, momentum_z = True)
    
    zz = (d1[:,0]-np.mean(d1[:,0]))*1e3
    II = d1[:,-2]; print(II.max())
    ee = d1[:,3]*1e6
    delgam = d1[:,2]*g_mec2*1e3

    fig, ax = plt.subplots()
    ax.plot(zz, II, 'r-')
    ax.plot([-1], [-1], 'b--')
    ax.set_xlabel(r'$\xi$ (mm)')
    ax.set_ylabel(r'Current (A)', color = 'r')
    #ax.set_xlim(-5, 5)
    #ax.set_ylim(0, 200)
    ax.tick_params(axis='y', colors='r')
    ax.grid()
    ax.legend(['current', 'emittance'])

    ax1 = ax.twinx()
    ax1.plot(zz, ee, 'b--')
    ax1.set_ylabel(u_emi_x, color = 'b')
    #ax1.set_ylim(0., 10)
    #ax1.set_yticks([0, 2, 4, 6, 8, 10])
    ax1.tick_params(axis='y', colors='b')

    fig.savefig('slice@'+fname+'.eps')
    
    
    fig, ax = plt.subplots()
    ax.plot(zz, II, 'r-')
    ax.plot([-1], [-1], 'b--')
    ax.set_xlabel(r'$\xi$ (mm)')
    ax.set_ylabel(r'Current (A)', color = 'r')
    #ax.set_xlim(-5, 5)
    #ax.set_ylim(0, 200)
    ax.tick_params(axis='y', colors='r')
    ax.grid()
    ax.legend(['current', 'emittance'])

    ax1 = ax.twinx()
    ax1.plot(zz, delgam, 'b--')
    ax1.set_ylabel(r'$\sigma_P$ (keV)', color = 'b')
    ax1.set_ylim(0., 50)
    #ax1.set_yticks([0, 2, 4, 6, 8, 10])
    ax1.tick_params(axis='y', colors='b')

    fig.savefig('slice2@'+fname+'.eps')
    
    #sys.exit()
    