
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 20:34:08 2022

@author: lixiangk
"""


from interface import *

# Get the directory of the currently running script
scriptdir = os.path.dirname(os.path.abspath(__file__))

# Add this directory to sys.path if it's not already there
if scriptdir not in sys.path:
    sys.path.insert(0, scriptdir)

# Now you can import the module as usual
from ObjOcelot import *


#%%% Re-scale the Twiss parameters
fname = 'beam_112A_2.0nC.ini'
fname = '368A.2809.002.1000'

betax, betay =  20.0,  0.75
alphax, alphay = 10.92*0.667, 3.25 

xrms, yrms = 1.5, 0.25
xx, yy = 1.0, 0.3

xemit, yemit = 8, 5
xemit, yemit = 5.5*1.5, 5.5*1.5

##
betax, betay =  7.934690,  0.4165
alphax, alphay = 5.257*0.65, 1.894


# alphax = 5

xrms, yrms = 1.25, 0.3
#xx, yy = 1.0, 0.3

xemit, yemit = 5, 5
###


P0 = 17
gamma = momentum2gamma(P0)

# betax = xrms**2/(xemit/gamma)
# betay = yrms**2/(yemit/gamma)

#alphax = -xx/(xemit/gamma)*1.45
#alphay = -yy/(yemit/gamma)

fout = 'matched@'+fname

kwargs = {}
kwargs.update(
    inputName = fname,
    outputName = fout,
    betax = betax,
    betay = betay,
    alphax = alphax,
    alphay = alphay,
    flipZ = False)
Matching(**kwargs)


workdir = r'/afs/ifh.de/group/pitz/data/lixiangk/sim3/2024/THzTransportDemo/_L-6.1ps-D-3.50mm-E1-58.19MV_m-phi1--10.28deg-E2-13.34MV_m-phi2--23.42deg-I-367.35A/'
os.chdir(workdir)
#%% High1.Scr1 to Undulator
Lquad, Rquad = 0.0675, 0.0215
rho, theta = 0.3, np.pi*60/180
astra_to_sco = 1.590444459638447

names = ['HIGH1.Q4', 'HIGH1.Q6', 'HIGH1.Q7'] 
grads = [-7.000000E-01,   9.675679E-01,  -4.564253E-01]

names += ['PST.QM3', 'PST.QT3', 'PST.QT6']
grads += [3.500000E-01,  -6.458153E-01,   3.426756E-01] # 002

names += ['HIGH2.Q3', 'HIGH2.Q4', 'HIGH2.Q5',
          'High3.Q1', 'High3.Q2', 'High3.Q3']
grads += [0.31912, -0.62146,  0.44067,
                  1.3814,  -3.04205,  1.72822] # matched, 002


grads = np.array(grads)/astra_to_sco

Distribution = 'ast_2.0nC.0528.003'
Distribution = 'ast.0528.001'
Zstart = 5.28
Zstop = 28.09

Run = 1

p_array_i = astraBeam2particleArray(filename=Distribution)
P0 = p_array_i.E*1e3

quadNames = names
nquads = len(quadNames)

quadGrads = grads
k = G2K(quadGrads, P0).T[:]

z0 = Zstart
cell = [Marker()]
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

#%
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

particleArray2astraBeam(p_array, filename="oce.%04.0f.%03d" % (Zstop*100, Run))
res = plot_twiss(tws_track, z0 = Zstart, label = 'transport-High1-to-HIGH3@%03d' % Run)

#%% In undulator

optics = 'undulator'

und = Undulator(3e-2, 113, 3.49)

l0 = (3.6-113*3e-2)/2

und_ent = pitz['High3.Und']
und_to_s2 = pitz['High3.Scr2']-und_ent-3.6
und_to_s2 = 0

flag = '003'
label = '%s@%s' % (optics, flag)

Distribution = 'oce_3.0nC.2809.001'
p_array_i = astraBeam2particleArray(filename=Distribution)

# Scale the charge
# p_array_i.q_array *= 1e-3

cell = [Marker(), Drift(l0), und, Drift(l0), Drift(und_to_s2), Marker()]
lat = MagneticLattice(cell, method=method)

navi = Navigator(lat)
# add physics processes from the first element to the last of the lattice
navi.add_physics_proc(sc1, lat.sequence[0], lat.sequence[-1])
# definiing of unit step in [m]
navi.unit_step = 0.02*1

# deep copy of the initial beam distribution 
p_array = deepcopy(p_array_i)
start = time.time()
tws_track, p_array = track(lat, p_array, navi); 
print(len(tws_track))
print("time exec: ", time.time() - start, "sec")

filename = 'oce_2.0nC.3167.%s' % (flag)
particleArray2astraBeam(p_array, filename=filename)
res = plot_twiss(tws_track, z0 = 0, label = label)
