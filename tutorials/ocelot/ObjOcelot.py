import sys
if sys.platform == "linux" or sys.platform == "linux2":
    sys.path.append(r'/afs/ifh.de/group/pitz/data/lixiangk/work/sync/python')
    rootdir = r'//afs/ifh.de/group/pitz/data/lixiangk/work'
elif sys.platform == "win32":
    sys.path.append(r'\\afs\ifh.de\group\pitz\data\lixiangk\work\sync\python')
    rootdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\work'
elif sys.platform == "darwin":
    print("OS X")
else:
    print("Unknown platform!")

#from PITZ_constants import *
#from Transport import *

import sys

import timeit
import pickle
import numpy as np

import os, shutil

# this python library provides generic shallow (copy) and deep copy (deepcopy) operations 
from copy import deepcopy
# import from Ocelot main modules and functions
from ocelot import *
# import from Ocelot graphical modules
from ocelot.gui.accelerator import *
# load beam distribution
# this function convert CSRtrack beam distribution to Ocelot format
# - ParticleArray. ParticleArray is designed for tracking.
# in order to work with converters we have to import 
# specific module from ocelot.adaptors
#from ocelot.adaptors.csrtrack2ocelot import *
from ocelot.adaptors.astra2ocelot import *

from interface import *

#astra_to_sco = 1.590444459638447
def twiss2rms(twiss, gamma = 18.954/0.511+1):
    emit_x = twiss[0]/gamma
    emit_y = twiss[1]/gamma
    
    rms_x = np.sqrt(twiss[2]*emit_x)
    rms_y = np.sqrt(twiss[3]*emit_y)
    cov_x = -twiss[4]*emit_x
    cov_y = -twiss[5]*emit_y
    
    rms = [rms_x, rms_y, cov_x, cov_y]
    return rms

def get_obj1(x, fname = 'BeamDynamics.dat', direc = None, **kwargs):

    if direc is not None:
        os.chdir(direc)
    
    diag = np.loadtxt(fname)
    # beta_x == beta_y
    obj =  [1-diag[-1,3]/diag[-1,4]]

    cmd = 'del BeamDynamics.dat'
    #os.system(cmd)

    if direc is not None:
        os.chdir('..')

    res = list(x)+list(obj)+list(diag[-1,:])
    filename = '.'+os.sep+'results_Q2_scan_Q1_%.2fT_m.dat' % x[0]
    with open(filename,'a') as f_handle:
        np.savetxt(f_handle,np.atleast_2d(res),fmt='%14.6E')

    return res

def get_obj2(x, fname = 'BeamDynamics.dat', direc = None, **kwargs):

    if direc is not None:
        os.chdir(direc)
    
    diag = np.loadtxt(fname)

    idx = -1
    #idx, _ = index_min(np.abs(diag[:,0]-diag[-1,0]+0.0675/2))
    print(len(diag), idx)
    
    # std_x == std_y
    obj =  [1-diag[idx,3]/diag[idx,4]]
    #obj =  [1-diag[-34,6]/diag[-34,7]]
    #obj =  [1-diag[-1,3]/diag[-1,4]]

    cmd = 'del BeamDynamics.dat'
    #os.system(cmd)

    if direc is not None:
        os.chdir('..')

    res = list(x)+list(obj)+list(diag[-1,:])
    filename = '.'+os.sep+'results_Q2_scan_Q1_%.2fT_m.dat' % x[0]
    with open(filename,'a') as f_handle:
        np.savetxt(f_handle,np.atleast_2d(res),fmt='%14.6E')

    return res

def get_obj3(x, fname = 'BeamDynamics.dat', direc = None, zmin = 0, **kwargs):

    if direc is not None:
        os.chdir(direc)
    
    diag = np.loadtxt(fname)

    if zmin == 0:
        n1 = -5
    else:
        n1 = len(diag)-np.sum(diag[:,0]>zmin)+1
        print(n1)
    
    # std_x == std_y
    tmp = diag[n1:,3]-diag[n1:,4]
    
    l0 = len(tmp)
    l1 = np.sum(tmp>=0)
    l2 = np.sum(tmp<=0)
    
    
    if l1 == l0 or l2 == 0:
        obj = [np.sum(tmp)/l0]
    elif tmp[1]>0:
        obj = [(np.sum(tmp[:l1])-np.sum(tmp[l1:]))/l0]
    else:
        obj = [(np.sum(tmp[:l2])-np.sum(tmp[l2:]))/l0]
        
    #tmp1 = np.sum(np.abs(diag[n1:,3]-diag[n1:,4]))
    #tmp2 = np.sum((diag[n1:,3]-diag[n1:,4]))
    #obj = [tmp2] # compare rms size
    #obj = [np.sum(diag[-100:,6]-diag[-100:,7])] # compare beta

    cmd = 'del BeamDynamics.dat'
    #os.system(cmd)

    if direc is not None:
        os.chdir('..')
    
    res = list(x)+list(obj)+list(diag[-1,:])
    filename = '.'+os.sep+'results_Q3_scan_Q1_%.2fT_m.dat' % x[0]
    with open(filename,'a') as f_handle:
        np.savetxt(f_handle,np.atleast_2d(res),fmt='%14.6E')
        
    return res

def ObjOcelotLattice(grads, positions = None, quadNames = None, Zstart = 0, 
                     Zstop = None, Lquad = 0.0675, P0 = 17.05, **kwargs):
    
    grads = np.array([grad/astra_to_sco for grad in grads])
    
    if positions is not None:
        nquads = len(positions)
    else:
        nquads = len(quadNames)
    
    quadGrads = grads
    k = G2K(quadGrads, P0).T[:]
    
    cell = [Marker()]
    
    z0 = Zstart
    for i in np.arange(nquads):
        
        if positions is not None:
            z1 = positions[i]
        else:
            z1 = pitz[quadNames[i]]
        
        if z1>=Zstop:
            break
        
        q1 = Quadrupole(l = Lquad, k1 = k[i])
        cell += [Drift(l = z1-z0-Lquad/2), q1]
        
        z0 = z1+Lquad/2
    
    if Zstop is None:
        Zstop = z0+1
        
    cell += [Drift(l = Zstop-z0), Marker()]
    
    return cell


def ObjOcelot(grads, positions = None, quadNames = None, Distribution = None,
              Zstart = 0, Zstop = None, step = 1, Lquad = 0.0675, save = False,
              obj_fun = get_obj2, Distribution1 = None, **kwargs):
    '''
    #quadNames = ['HIGH1.Q4', 'HIGH1.Q6', 'HIGH1.Q7'] 
    #grads = [-7.000000E-01,   9.675679E-01,  -4.564253E-01]
    #grads = np.array(grads)/astra_to_sco
    
    #Distribution = '368A.0528.002'
    #Zstart = 5.28
    #Zstop = 28.087#16.262
    #Run = 1
    
    Parameters
    ----------
    grads : TYPE
        DESCRIPTION.
    positions : TYPE, optional
        DESCRIPTION. The default is None.
    quadNames : TYPE, optional
        DESCRIPTION. The default is None.
    Distribution : TYPE, optional
        DESCRIPTION. The default is None.
    Lquad : TYPE, optional
        DESCRIPTION. The default is 0.0675.
    Zstop : TYPE, optional
        DESCRIPTION. The default is None.
    Zstart : TYPE, optional
        DESCRIPTION. The default is 0.
    beamline : TYPE, optional
        DESCRIPTION. The default is None.
    obj_fun : TYPE, optional
        DESCRIPTION. The default is get_obj2.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    obj : TYPE
        DESCRIPTION.

    '''
    
    x = grads
    grads = np.array([grad/astra_to_sco for grad in grads])
    
    print(Distribution)
    if Distribution is None:
        Distribution = 'ast.2550.001'
        
    # if Zstop is None:
    #     Zstop = positions[-1]+1
    # Zstop -= Zstart
    
    direc = ''
    for i in np.arange(len(x)):
        direc += str.format('Q%d-%.2fT_m-' %  (i+1, x[i]))
    direc = direc[:-1]; print(direc)

    CreateFolder(direc)
    #os.system(cmd)

    try:
        os.mkdir(direc)
    except IOError as err:
        print(err)
    
    os.chdir(direc)

    ######################
    p_array_i = astraBeam2particleArray(filename=Distribution)
    P0 = p_array_i.E*1e3
    
    cell = [Marker()]
    
    #quadNames = names
    if positions is not None:
        nquads = len(positions)
    else:
        nquads = len(quadNames)
    
    quadGrads = grads
    k = G2K(quadGrads, P0).T[:]
    
    z0 = Zstart
    for i in np.arange(nquads):
        
        if positions is not None:
            z1 = positions[i]
        else:
            z1 = pitz[quadNames[i]]
        
        if z1>Zstop:
            break
        
        q1 = Quadrupole(l = Lquad, k1 = k[i])
        cell += [Drift(l = z1-z0-Lquad/2), q1]
        
        z0 = z1+Lquad/2
        
    if Zstop is None:
        Zstop = z0+1
        
    cell += [Drift(l = Zstop-z0), Marker()]
    
    # initialization of tracking method
    method = MethodTM()
    # for second order tracking we have to choose SecondTM 
    method = {'global': SecondTM}
    # for first order tracking uncomment next line
    method = {'global': TransferMap}
    
    # lattice
    lat = MagneticLattice(cell, method=method)
    
    # space charge setting
    sc1 = SpaceCharge()
    sc1.nmesh_xyz = [31, 31, 31]
    sc1.step = step
    
    navi = Navigator(lat)
    # add physics processes from the first element to the last of the lattice
    navi.add_physics_proc(sc1, lat.sequence[0], lat.sequence[-1])
    # definiing of unit step in [m]
    navi.unit_step = 0.02*5
    
    # deep copy of the initial beam distribution 
    p_array = deepcopy(p_array_i)
    start = time.time()
    tws_track, p_array = track(lat, p_array, navi); 
    #print(len(tws_track))
    #print("time exec: ", time.time() - start, "sec")
    
    if save:
        particleArray2astraBeam(p_array, filename = Distribution1)
    res = plot_twiss(tws_track, z0 = Zstart, label = '', fout = 'BeamDynamics.dat')
    
    #######################
    
    os.chdir('..')
    
    obj = obj_fun(x, direc = direc, **kwargs)

    return obj

def plot_twiss(tws_track, z0 = 0, nrows = 1, label = '', fig_ext = '.png', fout = None, plot = True):
    
    bg = (tws_track[0].E*1e3)/0.51099895000; print(bg)
        
    res = []
    for i, tw in enumerate(tws_track):
        res.append([tw.s, tw.emit_x*bg*1e6, tw.emit_y*bg*1e6, 
                    tw.xx**0.5*1e3, tw.yy**0.5*1e3, tw.tautau**0.5*1e3,
                    tw.beta_x, tw.beta_y, tw.alpha_x, tw.alpha_y])
    res = np.array(res)
    res[:,0] += z0
    
    if fout == None:
        fout = 'track-%s.dat' % label
    np.savetxt(fout, res, fmt = '%14.6E')
    
    if plot:
        fig, ax = plt.subplots(nrows = nrows, ncols = 4//nrows, figsize = (3*4/nrows, 3*nrows))
        ax = ax.flatten()
        
        ax[0].plot(res[:,0], res[:,1], 'r-*')
        ax[0].plot(res[:,0], res[:,2], 'b-<')
        
        ax[1].plot(res[:,0], res[:,3], 'r-*')
        ax[1].plot(res[:,0], res[:,4], 'b-<')
        ax[1].plot(res[:,0], res[:,5], 'g-<')
        
        ax[2].plot(res[:,0], res[:,6], 'r-*')
        ax[2].plot(res[:,0], res[:,7], 'b-<')
        
        ax[3].plot(res[:,0], res[:,8], 'r-*')
        ax[3].plot(res[:,0], res[:,9], 'b-<')
        
        ylabels = ['Norm. emit. (um)', 'RMS size (mm)', 'beta function (m)', 'alpha function']
        titles = ['Emittance', 'RMS size', 'beta', 'alpha']
        for i in np.arange(len(ax)):
            ax[i].set_title(titles[i] + ' vs s')
            ax[i].set_xlabel(r'$s$ (m)')
            ax[i].set_ylabel(ylabels[i])
            #ax[i].grid()
        plt.tight_layout()
            
        fig.savefig(('emit-and-beta-function-along-beamline-%s' % label)+fig_ext)
    return res