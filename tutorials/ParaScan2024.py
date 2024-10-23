import sys
if sys.platform == "linux" or sys.platform == "linux2":
    sys.path.append(r'/afs/ifh.de/group/pitz/data/lixiangk/work/sync/python')
elif sys.platform == "win32":
    sys.path.append(r'\\afs\ifh.de\group\pitz\data\lixiangk\work\sync\python')
elif sys.platform == "darwin":
    print("OS X")
else:
    print("Unknown platform!")

#from _differentialevolution99 import *
#from my_object import *

from interface import *
from PITZ import *

from timeit import default_timer
import numpy as np

import time

from scipy.interpolate import interp1d, interp2d

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

# workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim3\2024\THzSim'
# os.chdir(workdir)

#%% Defind post-processing function
import difflib

def find_closest_folder(target_folder_name, search_directory = '.'):
    # Get a list of folders in the specified directory
    folders = [f for f in os.listdir(search_directory) if os.path.isdir(os.path.join(search_directory, f))]
    
    # Use difflib to find the closest match
    closest_match = difflib.get_close_matches(target_folder_name, folders, n=1)
    
    # Return the closest match or None if no match is found
    return closest_match[0] if closest_match else None

def post_ParaScan(x, *args):
    '''
    Photocathode to the EMSY1 at 5.277 m
    Parameters:
      x: an array or list of the variables to be optimized
    Returns:
      energy: the quantity that is envisoned as the "energy" of the sample
    '''
    
    BSA = x[1]
    
    sigma_x = sigma_y = BSA/4.
    phi_gun, phi_booster = x[3], x[5]
    Imain = x[6]
    
    MaxE_gun = x[2]
    phi_gun, phi_booster = x[3], x[5]
    #MaxE_gun = get_MaxE_gun(phi_gun, 6.3)
    
    MaxE_booster = get_MaxE_booster(MaxE_gun, phi_gun, phi_booster, 17)
    MaxB = I2B(Imain)
    
    Q_total = x[0]/1e3
    #Ipart = int(Q_total*50e3)

    direc = str.format('Q-%.2fpC-D-%.2fmm-E1-%.2fMV_m-phi1-%.2fdeg-E2-%.2fMV_m-phi2-%.2fdeg-I-%.2fA' %
                       (np.abs(x[0]), BSA, MaxE_gun, phi_gun, MaxE_booster, phi_booster, Imain))
        
    direc = find_closest_folder(direc);print(direc)
    print(direc)
    
    try:
        
        os.chdir(direc)
        try:
            fname = 'ast.0528.001'
            #fname = 'ast.2809.002'
            diag = BeamDiagnostics(fname = fname, energy = False)
        except:
            fname = 'ast.0527.001'
            diag = BeamDiagnostics(fname = fname, energy = False)
            
        res = list(x) + list(diag.x)
        os.chdir('..'); #os.system('rm -r '+direc)

        with open('ParaScan_1nC_3.5mm@5.28m.dat','a') as f_handle:
            np.savetxt(f_handle,np.atleast_2d(res),fmt='%14.6E')

    except:
        # if no output beam file, then set the objective to 10000
        os.chdir('..')
        print('Error: not accelerated well!')

    return 0

#x = [5, 0, 0, 350]
#obj = obj_fun(x)
#exit()

#%% Define goal function for Astra
def obj_ParaScan(x, *args):
    
    '''
    Created on Nov 17, 2019
    Simulation from photocathode via EMSY1 at 5.277 m to 20 m
    Goal function is a combination of average emittance and correlated energy spread exterpolated at undulator center, 29.25 m
    Laser spot fixed at 4 mm, gun phase at MMMG
    Variables to be optimized: gun and booster phases, solenoid current
    Parameters:
      x: an array or list of the variables to be optimized
    Returns:
      energy: the quantity that is envisoned as the "energy" of the sample
    '''
    
    
    Ipart = 250000
    
    Q_total = -x[0]/1e3 # to nC
    
    BSA = x[1]
    
    FWHM = 7e-3 # for Gaussian
    #FWHM = 20e-3 # for Flattop
    #FWHM = 250e-6 # for Pharos shortest pulse
    #FWHM = x[0]*1e-3
    
    # uniform distribution
    sigma_x = sigma_y = BSA/4.
    C_sigma_x = C_sigma_y = 0
    
    #sigma_x = sigma_y = 3.0/2.355 # Gaussian truncated
    #C_sigma_x = C_sigma_y = BSA/2./sigma_x
    
    sigma_x = sigma_y = 0.83 #0.96 # Gaussian truncated
    C_sigma_x = C_sigma_y = BSA/2./sigma_x
    
    phi_gun, phi_booster = x[3], x[5]
    #MaxE_gun = get_MaxE_gun(phi_gun, 6.3)
    #MaxE_gun = x[2]
    
    MaxE_booster = get_MaxE_booster(MaxE_gun, phi_gun, phi_booster, 17)
    #MaxE_booster = x[4]
    
    Imain = x[6]
    MaxB = I2B(Imain)
    
    field_maps = rootdir+os.sep+'sync'+os.sep+'field-maps'
    
    #Distribution = 'beam.ini'
    Distribution = f'BSA-{BSA:.1f}mm.ini'
    #Distribution = 'ast.0528.011'
    generator = Generator1(FNAME = Distribution, IPart = Ipart, Species = 'electrons', Q_total = Q_total,
                           Ref_Ekin = 0.0e-6, LE = 0.55e-3*1, dist_pz = 'i',
                           Dist_z = 'g', sig_clock = FWHM/2.355, Cathode = True,
                           Dist_x = '2', sig_x = sigma_x, Dist_px = 'g', Nemit_x = 0,
                           Dist_y = '2', sig_y = sigma_y, Dist_py = 'g', Nemit_y = 0,
                           C_sig_x = C_sigma_x, C_sig_y = C_sigma_y) # Gaussian by default
    
    #generator.set(Ref_Ekin = 1e-6, LE = 0, dist_pz = 'g', Nemit_x = 0.17, Nemit_y = 0.18)
    #generator.set(Dist_x = 'r', sig_x = BSA/4.0, Dist_y = 'r', sig_y = BSA/4.0) # Transverse uniform
    #generator.set(Dist_z = 'p', Lt = FWHM, Rt = 1e-3) # Temporal flattop
    
    Run = 1
    newrun = Module('Newrun', Run = Run, Head = 'PITZ beam line simulation', Distribution = Distribution, CathodeS = True,
                    Auto_Phase = True, Track_All = True, check_ref_part = False, Lprompt = False, Max_step=200000)
    #newrun.set(Xoff = 0.5, Yoff = 0.5)
    #newrun.set(Qbunch = Q_total)
    #newrun.set(Run = 1, XYrms = sigma_x)
    
    charge = Module('Charge', LSPCH = True, Lmirror = True, Nrad = 40, Nlong_in = 50, N_min = 50, Max_scale = 0.05, Max_count = 20)
    charge.set(N_min= 50)
    #charge.set(L2D_3D = True, Z_TRANS = 5, NXF = 32, NYF = 32, NZF = 32, min_grid_trans = 0.03e-3)
    #charge.set(LSPCH3D = True, NXF = 32, NYF = 32, NZF = 32)
    #charge.set(LSPCH = False)
    
    cavity = Module('Cavity', LEfield = True, File_Efield = 
                    [field_maps+os.sep+gun_profile, field_maps+os.sep+'CDS14_15mm.txt'],
                    MaxE = [MaxE_gun, MaxE_booster], C_pos = [0., 2.675], Nue = [1.3, 1.3], Phi = [phi_gun, phi_booster])
    
    soleno = Module('Solenoid', LBfield = True, File_Bfield = [field_maps+os.sep+'gunsolenoidsPITZ.txt'],
                    MaxB = [-MaxB], S_pos = [0.], S_xrot = [0*1e-3], S_yrot = [0])
    
    Zstart, Zstop = 0, 5.28
    Zemit = int((Zstop-Zstart)/0.01)
    output = Module('Output', Zstart = Zstart, Zstop = Zstop, Zemit = Zemit, Zphase = 1, RefS = True, EmitS = True,
                    PhaseS = True, TrackS = False, LandFS = True, C_EmitS = True, LPROJECT_EMIT = True,
                    LOCAL_EMIT = False, Screen = [5.28])
    #output.set(Zstop = 0.5, Zemit = 50) # emission
    
    apertu = Module('Aperture', LApert = True, File_Aperture = [field_maps+os.sep+'app.txt'])
    
    #Qsol = 0 #-0.0034*x[-2]/365.0
    #Qsol_zrot = 0 # -0.4363 #x[-1]/180.0*np.pi
    #quadru = Module('Quadrupole', LQuad = True, Q_type = ['../Qsol.data'], Q_pos = [0], Q_grad = [Qsol], Q_zrot = [Qsol_zrot])
    
    astra = Astra()
    astra.add_modules([newrun, charge, cavity, soleno, output])
    
    direc = str.format('Q-%.2fpC-D-%.2fmm-E1-%.2fMV_m-phi1-%.2fdeg-E2-%.2fMV_m-phi2-%.2fdeg-I-%.2fA' %\
                       (np.abs(x[0]), BSA, MaxE_gun, phi_gun, MaxE_booster, phi_booster, Imain))
    
    #direc = str.format('D-%.2fmm-I-%.0fA' % (x[1], x[-1]))
    ###
    os.system('mkdir -p '+direc)
    #os.chdir(direc)
    
    job_name = direc+'-%03d' % Run
    #job_name = 'myjob.'+('BSA-%.2fmm-I-%.2fA' % (x[0], x[3]))
    gen_name = 'gen' #+`Run`
    ast_name = 'ast' #+`Run`
    
    generator.write(direc+os.sep+gen_name+'.in')
    astra.write(direc+os.sep+ast_name+'.in')
    #astra.qsub(job_name, ast_name, gen_name, direc)
    #os.system('qsub '+job_name+'.sh')
    
    astra.submit(job_name, ast_name, gen_name, direc)
    
    return
    ###
    
    os.system('mkdir -p '+direc)
    os.chdir(direc)
    
    generator.write('gen.in')
    astra.write('ast.in')
    
    os.system('generator gen.in 2>&1 | tee gen.log')
    os.system('astra ast.in 2>&1 | tee ast.log')
    
    return

#%% Job submission and cancaling

## for condor jobs
# alias "cq=condor_q"
# alias to_sub_all="ls *.submit         | awk '{print \"condor_submit \"\$0}'"
# alias to_del_all="condor_q | grep \"ID:\" | awk '{print(\"condor_rm \", \$3)}'"
# to_sub_all | bash -

#%% Parameter scan using Astra

# Charge, pC
var0 = [1000]

# BSA, mm
var1 = [3.5]

# Gun phase, degree
var2 = np.linspace(-9, 0, 4)
var2 = [0]

# Booster phase, degree
var3 = np.linspace(-36, -16, 6)
#var3 = [-28, -24, -20]

# Solenoid, A
var4 = np.linspace(360, 374, 8)

combi = np.array([[v0, v1, 57.55, v2, 14.5, v3, v4] for v0 in var0 for v1 in var1 for v2 in var2 for v3 in var3 for v4 in var4])

for x in combi:
    print(x)
    #obj_ParaScan(x)
    post_ParaScan(x)
    pass

#%% Stop here if launched from terminal
exit()

#%% Parameter scan using Ocelot
def G2K(G, P0, qn = 1.):
    '''
    Parameters
      G: gradient, T/m
      P0: momentum, MeV/c
      qn = q/qe: number of charge
    Returns
      K: focal strength, m^-2
    '''
    return G*qn/(P0*1e6/g_c)

#%%% Either define a function for parameter scan
try:
    def obj_OcelotScan(x, *args, **kwargs):
        '''
        Define optics and run simulation with ocelot.
        The map of name and position pairs of magnets (mainly quadrupoles) is 
        neccessary.
    
        Parameters
        ----------
        x : 1D array or list 
            Gradients of quadrupoles, T/m.
        *args : 1D list
            Names of quadrupoles, e.g., High1.Q1, High1.Q2.
        **kwargs : 
            Optional input arguments: 
                Zstart, Zstop, Distribution, Run, Step ...
                mostly the naming convention of Astra is considered
    
        Returns
        -------
        res: 1D array or list
            The beam parameters (statistics) in the end: [z, eps_x^n, eps_y^n,
            sigma_x, sigma_y, beta_x, beta_y, alpha_x, alpha_y]
        tws_track, p_array:
            Twiss parameters and distribution of last step.
    
        '''
        
        
        # Initiate the arguments that are acceptable from kwargs
        Zstart = 0
        Zstop = pitz['High1.Scr1']
        
        Distribution = 'beam.ini'
        Run = 1
        
        Step = 1 # Step for space charge calculation in Ocelot
        Direc = '.'
        
        if len(kwargs)>0:
            if 'Zstart' in list(kwargs.keys()):
                Zstart = kwargs['Zstart']
    
            if 'Distribution' in list(kwargs.keys()):
                Distribution = kwargs['Distribution']
            
            if 'Zstop' in list(kwargs.keys()):
                Zstop = kwargs['Zstop']
    
            if 'Run' in list(kwargs.keys()):
                Run = kwargs['Run']
        
            if 'Step' in list(kwargs.keys()):
                Step = kwargs['Step']
                
            if 'Direc' in list(kwargs.keys()):
                Direc = kwargs['Direc']
        
        cwd = os.getcwd();
        os.chdir(Direc); print(Direc)
        
        p_array_i = astraBeam2particleArray(filename=Distribution)
        P0 = p_array_i.E*1e3 # to MeV/c
        
        cell = [Marker()]
        if len(args)>0:
            quadNames = args[0]
            nquads = len(quadNames)
            
            if len(args)>1:
                quadGrads = args[1]
            else:
                quadGrads = np.zeros(nquads)
            k = G2K(quadGrads, P0).T[:]
            
            z0 = Zstart
            for i in np.arange(nquads):
                
                z1 = pitz[quadNames[i]]
                q1 = Quadrupole(l = Lquad, k1 = k[i])
                cell += [Drift(l = z1-z0-Lquad/2), q1]
                
                z0 = z1
            cell += [Drift(l = Zstop-z0-Lquad/2), Marker()]
            
        else:
            cell += [Drift(l = Zstop-Zstart), Marker()]
            
        # initialization of tracking method
        method = MethodTM()
        # for second order tracking we have to choose SecondTM 
        method.global_method = SecondTM
        # for first order tracking uncomment next line
        # method.global_method = TransferMap
    
        # lattice
        lat = MagneticLattice(cell, method=method)
    
        # space charge setting
        sc1 = SpaceCharge()
        sc1.nmesh_xyz = [31, 31, 31]
        sc1.step = Step
    
        navi = Navigator(lat)
        # add physics processes from the first element to the last of the lattice
        navi.add_physics_proc(sc1, lat.sequence[0], lat.sequence[-1])
        # definiing of unit step in [m]
        navi.unit_step = 0.02
    
        # deep copy of the initial beam distribution 
        p_array = deepcopy(p_array_i)
        start = time.time()
        tws_track, p_array = track(lat, p_array, navi); 
        #print(len(tws_track))
        #print("time exec: ", time.time() - start, "sec")
        
        filename = f"oce.{Zstop*100:04.0f}.{Run:03d}"
        particleArray2astraBeam(p_array, filename = filename)
        try:
            res = plot_twiss(tws_track, z0 = Zstart, label = 'transport@{Run:03d}', plot = True)
        except Exception as err:
            print(err)
            print('Function `plot_twiss` may no be defined !')
            res = [None]
            
        os.chdir(cwd)
        
        return [res[-1], tws_track, p_array] 
except Exception as err:
    print(err)
#%%% Or run directly
names = ['HIGH1.Q4', 'HIGH1.Q6', 'HIGH1.Q7'] 
grads = [-7.000000E-01,   9.675679E-01,  -4.564253E-01]

names += ['PST.QM3', 'PST.QT3', 'PST.QT6']
grads += [3.500000E-01,  -6.458153E-01,   3.426756E-01]

names += ['HIGH2.Q3', 'HIGH2.Q4', 'HIGH2.Q5',
          'High3.Q1', 'High3.Q2', 'High3.Q3']
grads += [0.31912, -0.62146,  0.44067,
                  1.3814,  -3.04205,  1.72822] 

astra_to_sco = 1.590444459638447
grads = np.array(grads)/astra_to_sco

Distribution = 'ast.0528.001'
#Distribution = 'ast.0528.002' # 1 nC
Zstart = 5.28
Zstop = 28.09

Run = 1

p_array_i = astraBeam2particleArray(filename=Distribution)
P0 = p_array_i.E*1e3

cell = [Marker()]

quadNames = names
nquads = len(quadNames)

quadGrads = grads
k = G2K(quadGrads, P0).T[:]

z0 = Zstart
for i in np.arange(nquads):
    
    z1 = pitz[quadNames[i]]
    q1 = Quadrupole(l = Lquad, k1 = k[i])
    cell += [Drift(l = z1-z0-Lquad/2), q1]
    
    z0 = z1+Lquad/2
cell += [Drift(l = Zstop-z0), Marker()]

# initialization of tracking method
method = MethodTM()
# for second order tracking we have to choose SecondTM 
method.global_method = SecondTM
# for first order tracking uncomment next line
# method.global_method = TransferMap

# lattice
lat = MagneticLattice(cell, method=method)

# space charge setting
sc1 = SpaceCharge()
sc1.nmesh_xyz = [31, 31, 31]
sc1.step = 1

navi = Navigator(lat)
# add physics processes from the first element to the last of the lattice
navi.add_physics_proc(sc1, lat.sequence[0], lat.sequence[-1])
# definiing of unit step in [m]
navi.unit_step = 0.02*1

# deep copy of the initial beam distribution 
p_array = deepcopy(p_array_i)
start = time.time()
tws_track, p_array = track(lat, p_array, navi); 
#print(len(tws_track))
#print("time exec: ", time.time() - start, "sec")

particleArray2astraBeam(p_array, filename=f"oce.{Zstop*100:04.0f}.{Run:03d}")
res = plot_twiss(tws_track, z0 = Zstart, label = 'transport@')

#%%%    
# Charge, pC
var0 = [1000]

# BSA, mm
var1 = [3.5]

# Gun phase, degree
var2 = np.linspace(-9, 0, 4)
var2 = [0]

# Booster phase, degree
var3 = np.linspace(-36, -16, 6)
var3 = [-28, -24, -20]

# Solenoid, A
var4 = np.linspace(340, 372, 5)
var4 = np.linspace(360, 374, 8)

combi = np.array([[v0, v1, 57.55, v2, 14.5, v3, v4] for v0 in var0 
                  for v1 in var1 for v2 in var2 for v3 in var3 for v4 in var4])

for x in combi:
    
    # direc = 
    # runcell('Define macro script', '//afs/ifh.de/group/pitz/data/lixiangk/sim3/2024/THzSim/ParaScan2024.py')
    
    # Create the name of the directory
    direc = '.'
    
    names = ['HIGH1.Q4', 'HIGH1.Q6', 'HIGH1.Q7'] 
    grads = [-7.000000E-01,   9.675679E-01,  -4.564253E-01]

    names += ['PST.QM3', 'PST.QT3', 'PST.QT6']
    grads += [3.500000E-01,  -6.458153E-01,   3.426756E-01]

    names += ['HIGH2.Q3', 'HIGH2.Q4', 'HIGH2.Q5',
              'High3.Q1', 'High3.Q2', 'High3.Q3']
    grads += [0.31912, -0.62146,  0.44067,
                      1.3814,  -3.04205,  1.72822] 

    astra_to_sco = 1.590444459638447
    grads = np.array(grads)/astra_to_sco

    Distribution = 'ast.0528.001'
    #Distribution = 'ast.0528.002' # 1 nC
    Zstart = 5.28
    Zstop = 28.09

    Run = 1
    
    obj_OcelotScan(x, names, grads, Distribution = Distribution, Zstart = Zstart,
                  Zstop = Zstop, Run = Run, Direc = direc)
    pass

#%% Make plots after post-processing, first load data
import pandas as pd

data = np.loadtxt('ParaScan_2nC_3.5mm@5.28m.dat')
scanned = ['charge', 'BSA', 'E1', 'phi1', 'E2', 'phi2', 'Imain']

diaged = list(BeamDiagnostics().keyIndex.keys())
df = pd.DataFrame(data, columns = scanned+diaged)

print(df.keys())

flag = 'BSA'; var = [3.5]

#%%% Higher-order-energy-spread-vs-booster-phase-vs-gun-phase
flag= 'Imain'
select = np.abs(df[flag]-370)<1; print(select)
ds = df[select]

var = np.linspace(-6, 0, 3)
flag = 'phi1'

fig, ax = plt.subplots()
for v0 in var:
    select = np.abs(ds[flag]-v0)<0.01; print(select)
    ds1 = ds[select]
    # ax.plot(ds['Imain'], ds['nemit_x'], 'r')
    # ax.plot(ds['Imain'], ds['nemit_y'], 'b')
    ax.plot(ds1['phi2'], ds1['std_Ek2'])

ax.grid()
ax.set_xlabel(r'$\phi_{\rm boo}$ (degree)')
ax.set_ylabel(r'Higher order $\Delta P$ (keV/c)')

ax.legend(['%s = %.1f deg' % (flag, v0) for v0 in var])
#ax.legend(['%s = %.0f deg' % (flag, v0) for v0 in var])

fig.savefig('std_Ek2-vs-phi2-vs-phi1-370A@2nC@5.28m.png')

#%%% Energy-spread-vs-booster-phase-vs-gun-phase
flag= 'Imain'
select = np.abs(df[flag]-364)<1; print(select)
ds = df[select]

var = np.linspace(-6, 0, 3)
flag = 'phi1'

fig, ax = plt.subplots()
for v0 in var:
    select = np.abs(ds[flag]-v0)<0.01; print(select)
    ds1 = ds[select]
    # ax.plot(ds['Imain'], ds['nemit_x'], 'r')
    # ax.plot(ds['Imain'], ds['nemit_y'], 'b')
    ax.plot(ds1['phi2'], ds1['std_Ekin'])

ax.grid()
ax.set_xlabel(r'$\phi_{\rm boo}$ (degree)')
ax.set_ylabel(r'$\Delta P$ (keV/c)')

ax.legend(['%s = %.1f deg' % (flag, v0) for v0 in var])
#ax.legend(['%s = %.0f deg' % (flag, v0) for v0 in var])

fig.savefig('std_Ek-vs-phi2-vs-phi1-370A@2nC@5.28m.png')

#%%% Ipeak-vs-booster-phase-vs-gun-phase
flag= 'Imain'
select = np.abs(df[flag]-364)<1; print(select)
ds = df[select]

var = np.linspace(-6, 0, 3)
flag = 'phi1'

fig, ax = plt.subplots()
for v0 in var:
    select = np.abs(ds[flag]-v0)<0.01; print(select)
    ds1 = ds[select]
    # ax.plot(ds['Imain'], ds['nemit_x'], 'r')
    # ax.plot(ds['Imain'], ds['nemit_y'], 'b')
    ax.plot(ds1['phi2'], ds1['I1'])

ax.grid()
ax.set_xlabel(r'$\phi_{\rm boo}$ (degree)')
ax.set_ylabel(r'Ipeak (A)')

ax.legend(['%s = %.1f deg' % (flag, v0) for v0 in var])
#ax.legend(['%s = %.0f deg' % (flag, v0) for v0 in var])

fig.savefig('Ipeak-vs-phi2-vs-phi1-370A@2nC@5.28m.png')

#%%% Xrms-vs-booster-phase-vs-gun-phase
flag= 'Imain'
select = np.abs(df[flag]-364)<1; print(select)
ds = df[select]

var = np.linspace(-9, 0, 4)
flag = 'phi1'

fig, ax = plt.subplots()
for v0 in var:
    select = np.abs(ds[flag]-v0)<0.01; print(select)
    ds1 = ds[select]
    # ax.plot(ds['Imain'], ds['nemit_x'], 'r')
    # ax.plot(ds['Imain'], ds['nemit_y'], 'b')
    ax.plot(ds1['phi2'], ds1['std_x'])

ax.grid()
ax.set_xlabel(r'$\phi_{\rm boo}$ (degree)')
ax.set_ylabel(r'$\sigma_x$ (mm)')

ax.legend(['%s = %.1f deg' % (flag, v0) for v0 in var])
#ax.legend(['%s = %.0f deg' % (flag, v0) for v0 in var])

fig.savefig('std_x-vs-phi2-vs-phi1-364A@2nC@5.28m.png')

#%%% Qb-vs-charge-vs-BSA
var = [4]
flag = 'BSA'

fig, ax = plt.subplots()
for v0 in var:
    select = np.abs(1-df[flag]/v0)<0.05
    ds = df[select]
    #ds =df
    # ax.plot(ds['Imain'], ds['nemit_x'], 'r')
    # ax.plot(ds['Imain'], ds['nemit_y'], 'b')
    ax.plot(ds['charge'], -ds['Q_b']*1e12)

ax.grid()
ax.set_xlabel(r'BSA (mm)')
ax.set_xlabel(r'Imain (A)')
ax.set_ylabel(r'Bunch charge (pC)')

ax.legend(['%s = %.1f mm' % (flag, v0) for v0 in var])

fig.savefig('Qb-vs-charge-vs-BSA.png')

#%%% Qb-vs-Imain-vs-
fig, ax = plt.subplots()
for v0 in var:
    select = np.abs(1-df[flag]/v0)<0.05
    ds = df[select]
    # ax.plot(ds['Imain'], ds['nemit_x'], 'r')
    # ax.plot(ds['Imain'], ds['nemit_y'], 'b')
    ax.plot(ds['Imain'], -ds['Q_b']*1e12)

ax.grid()
ax.set_xlabel(r'Imain (A)')
ax.set_ylabel(r'Bunch charge (pC)')

ax.legend(['%s = %.1f deg' % (flag, v0) for v0 in var])

fig.savefig('Qb-vs-Imain-vs-BSA.png')

#%%% nemit-vs-Imain-vs-
s0 = np.abs(0-df['phi1'])<1

var = [-24, -22, -20]
var = [-20]
flag= 'phi2'
fig, ax = plt.subplots()
for v0 in var:
    select = np.abs(df[flag]-v0)<0.01; print(select)
    ds = df[select*s0]
    # ax.plot(ds['Imain'], ds['nemit_x'], 'r')
    # ax.plot(ds['Imain'], ds['nemit_y'], 'b')
    ax.plot(ds['Imain'], np.sqrt(ds['nemit_x']*ds['nemit_y']))

ax.grid()
ax.set_xlabel(r'Imain (A)')
ax.set_ylabel(r'Norm. emittance ($\mu$m)')

ax.legend(['%s = %.1f deg' % (flag, v0) for v0 in var])
#ax.legend(['%s = %.0f deg' % (flag, v0) for v0 in var])

#fig.savefig('nemit-vs-Imain-vs-BSA@2nC@3.5mm.png')

#%%% nemit-vs-Imain-vs-BSA
#select = np.abs(df['phi2']+5)<0.01; print(select)
#ds0 = df[select]

flag = 'BSA'
#var = [0.6, 0.8, 1.0, 1.2]
fig, ax = plt.subplots()
for v0 in var:
    select = np.abs(df[flag]-v0)<0.1; print(select)
    ds = df[select]
    # ax.plot(ds['Imain'], ds['nemit_x'], 'r')
    # ax.plot(ds['Imain'], ds['nemit_y'], 'b')
    ax.plot(ds['Imain'], np.sqrt(ds['nemit_x']*ds['nemit_y']))

ax.grid()
ax.set_xlabel(r'Imain (A)')
ax.set_ylabel(r'Norm. emittance ($\mu$m)')

ax.set_xlim(356, 376)
#ax.set_ylim(0, 3)

ax.legend(['%s = %.1f mm' % (flag, v0) for v0 in var])
#ax.legend(['%s = %.0f deg' % (flag, v0) for v0 in var])

fig.savefig('nemit-vs-Imain-vs-BSA@2nC.png')

#%%% rms-vs-Imain-vs-
s0 = np.abs(0-df['phi1'])<1

var = [-24, -22, -20]
var = [-20]
flag= 'phi2'
fig, ax = plt.subplots()
for v0 in var:
    select = np.abs(1-df[flag]/v0)<0.05
    ds = df[s0*select]
    # ax.plot(ds['Imain'], ds['nemit_x'], 'r')
    # ax.plot(ds['Imain'], ds['nemit_y'], 'b')
    ax.plot(ds['Imain'], np.sqrt(ds['std_x']*ds['std_y']))

ax.grid()
ax.set_xlabel(r'Imain (A)')
ax.set_ylabel(r'RMS size (mm)')

#ax.set_xlim(356, 376)
#ax.set_ylim(0, 3)

ax.legend(['%s = %.1f mm' % (flag, v0) for v0 in var])

#fig.savefig('rms-vs-Imain-vs-BSA@2nC@3.5mm.png')

#%%% cor_Ek-vs-phi2-vs-
#flag = 'Imain'
#v0 = 388
#select = np.abs(1-df[flag]/v0)<0.001; ds = df[select]

flag = 'BSA'; var = [1.8, 2]

fig, ax = plt.subplots()
for v0 in var:
    select = np.abs(1-df[flag]/v0)<0.01
    ds = df[select]
    #ds = df
    # ax.plot(ds['Imain'], ds['nemit_x'], 'r')
    # ax.plot(ds['Imain'], ds['nemit_y'], 'b')
    ax.plot(ds['phi2'], ds['cor_Ekin'])
    ax.plot(ds['phi2'], ds['std_Ekin'])

ax.grid()
ax.set_xlabel(r'$\phi_{boo}$ (degree)')
ax.set_ylabel(r'Energy spread (keV)')

ax.legend([r'Correlated', 'RMS'])

fig.savefig('cor_Ekin-vs-phi2-BSA.png')
#exit()

