# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 15:15:22 2023

@author: lixiangk
"""
from interface import *

def Halton(i, j):
    '''
    Hamseley sequence reproduced from Mithra2.

    Parameters
    ----------
    i : TYPE
        DESCRIPTION.
    j : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    if (i > 240):
        i = i%240
        print("Dimension can not be larger than 240. Residue is taken instead: %d" % i)
        
    prime = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
   
    prime = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 
             31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
             73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 
             127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 
             179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 
             233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 
             283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 
             353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 
             419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 
             467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 
             547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 
             607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 
             661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 
             739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 
             811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 
             877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 
             947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 
             1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 
             1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 
             1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 
             1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 
             1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 
             1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 
             1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511] # 240
    #int p0, p, k, k0, a;
    
    x = 0.0;
    k0 = j + 1;
    p = prime[i];

    p0 = p;
    k  = k0;
    x  = 0.0;
    while (k > 0):
        a = k%p;
        x += a/float(p0);
        k = int(k/p);
        p0 *= p;
      
    return 1.0 - x;
halton = Halton;

# def HaltonNorm(j, mu = 0, sig = 1):
#     r1 = Halton(8, j)
#     r2 = Halton(9, j)   
#     return sig*np.sqrt(-2*np.log(r1))*np.cos(2*np.pi*r2)+mu

def HaltonNorm(j, mu = 0, sig = 1, c_sig = 0):
    r1 = Halton(8, j)
    if c_sig == 0:
        rr = sp.special.erfinv(2*r1-1)*np.sqrt(2)*sig+mu
    else:
        ccc = sp.special.erf(c_sig/np.sqrt(2))
        rr = sp.special.erfinv(2*r1*ccc-ccc)*np.sqrt(2)*sig+mu
    return rr


erf = sp.special.erf
erfinv = sp.special.erfinv

def HaltonNorm2(j, mu = 0, sig = 1, a = -10, b = 10, i = 8):
    # F(x, a) = (erf(x/np.sqrt(2))-erf(a/np.sqrt(2)))/2/ccc -> 0, 1
    Fa = erf(a/np.sqrt(2))
    Fb = erf(b/np.sqrt(2))
    r1 = Halton(i, j)
    r2 = erfinv(2*r1*(Fb-Fa)/2+Fa)*np.sqrt(2)
    return r2

#
from scipy.stats import qmc
def HaltonN(n, iprime = 1):
    sampler = qmc.Halton(d = iprime, scramble = True)
    rn = sampler.random(n = n).reshape((n,iprime)); #print(rn.shape)
    return rn[:,iprime-1]

def HaltonNNorm(n, mu = 0, sig = 1, a = -10, b = 10, iprime = 8):
    # F(x, a) = (erf(x/np.sqrt(2))-erf(a/np.sqrt(2)))/2/ccc -> 0, 1
    Fa = erf(a/np.sqrt(2))
    Fb = erf(b/np.sqrt(2))
    
    r1 = HaltonN(n, iprime)
    r2 = erfinv(2*r1*(Fb-Fa)/2+Fa)*np.sqrt(2)
        
    return r2

def HaltonNNorm2(samples, mu = 0, sig = 1, a = -10, b = 10):
    # F(x, a) = (erf(x/np.sqrt(2))-erf(a/np.sqrt(2)))/2/ccc -> 0, 1
    Fa = erf(a/np.sqrt(2))
    Fb = erf(b/np.sqrt(2))
    
    #r1 = HaltonN(n, iprime)
    r1 = samples
    r2 = erfinv(2*r1*(Fb-Fa)/2+Fa)*np.sqrt(2)
        
    return r2

#%%
prime = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 
        31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
        73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 
        127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 
        179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 
        233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 
        283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 
        353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 
        419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 
        467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 
        547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 
        607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 
        661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 
        739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 
        811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 
        877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 
        947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 
        1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 
        1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 
        1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 
        1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 
        1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 
        1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 
        1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511]
from scipy.stats import qmc
sampler = qmc.Halton(d=1, scramble=False, seed = 100)
rr = sampler.random(n=100)


fig, ax = plt.subplots(figsize = (5, 4))
ax.hist(rr, bins = 20, histtype = r'step')

#rr = np.array([HaltonNorm2(j, a = -3, b = 0) for j in np.arange(10000)])
rr = np.array([Halton(100, j) for j in np.arange(100)])
#fig, ax = plt.subplots(figsize = (5, 4))
ax.hist(rr, bins = 20, histtype = r'step')
#%% Genesis1.3-version4 modify theta from existing output particles
import h5py

workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2021\Genesis-demo\Minimal'
workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2022\THzComm2\2nC\SASE\002-G4'
os.chdir(workdir)

fname = 'pitz.thz.1.0.par.h5'
fname = 'temp.par.h5'
f = h5py.File(fname, 'r+')

Qb = 2e-9
sigma_z = 2e-3
sigma_t = sigma_z/g_c
Ipeak = Qb/np.sqrt(2*np.pi)/sigma_t

lambda0 = 100e-6
k0 = 2*np.pi/lambda0

fG = lambda x, mu = 0, sig = sigma_z:Ipeak*np.exp(-(x-mu)**2/2/sig**2)

c_sig = 2
Ns_b = int(sigma_z*c_sig*2/lambda0)//2*2+2

#fig, ax = plt.subplots()

Ns = f['slicecount'][0]

Nm_slices = 0
stat = []
for islice in np.arange(0, Ns_b+20):
    name = str.format('slice%06d' % (islice+2))
    
    current = f[name+'/current'][0]
    if current > 0:
        theta = f[name+'/theta'][:]; Nm_slice = len(theta)
    # else:
    #     theta = f[name+'/theta'][:]*0; Nm_slice = len(theta)
        
        x1 = -c_sig*sigma_z+islice*lambda0
        a = (x1-0.5*lambda0)/sigma_z
        b = (x1+0.5*lambda0)/sigma_z
        
        rr = np.array([HaltonNorm2(Nm_slices+j, a = a, b = b) for j in np.arange(Nm_slice)])-a
        
        current1 = fG(x1)
        theta1 = 2*np.pi*(rr*sigma_z/lambda0)
        
        
        Nm_slices += Nm_slice
        
        print(islice, current, current1)
        
        bf2 = np.abs(np.mean(np.exp(1j*theta)))**2
        bf2_new = np.abs(np.mean(np.exp(1j*theta1)))**2
        stat.append([islice, x1, current, current1, bf2, bf2_new])
    
        f[name+'/current'][0] = current1
        del f[name+'/theta']
        dset = f.create_dataset(name+'/theta', data = theta1)

    #ax.hist(theta/k0+x1, bins = 200, histtype = r'step')
    #ax.hist(theta1/k0+x1, bins = 200, histtype = r'step')
f.close() 
    
stat = np.array(stat)

fig, ax = plt.subplots(figsize = (4, 6), nrows = 2, sharex = True)
ax1, ax2 = ax.flatten()
ax1.plot(stat[:,0], stat[:,2])
ax1.plot(stat[:,0], stat[:,3])

ax2.plot(stat[:,0], stat[:,4])
ax2.plot(stat[:,0], stat[:,5])

ax2.set_yscale('log')

#%%
workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2022\THzComm2\2nC\SASE\002-G4'
os.chdir(workdir)

fname = 'pitz.thz.2.out.h5'

pg = PostGenesis(fname, debug = 0, fig_ext = '.png')

zplot = pg.zplot
zpower = pg.zpower
zenergy = pg.zenergy

plt.figure(999)
plt.plot(zplot, zenergy*1e6)
plt.yscale('log')

workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2022\THzComm2\2nC\SASE\002-bunch'
os.chdir(workdir)

ext = ''
EE = np.loadtxt('energy-z-ipseed'+ext+'.txt'); EE = EE.T; print (EE.shape)

plt.plot(EE[:,0], EE[:,1]*1e6)







#%% Postprocess output
fname = 'pithz.1.out'
#fname = 'template.out'
pg = PostGenesis(fname, debug = 1, fig_ext = '.png')

#%% Postprocess output, batch
fig, ax = plt.subplots(figsize = (5, 4))
for i in [2, 3, 4, 5]:
    fname = 'pithz.%d.out' % i

    pg = PostGenesis(fname, debug = 0)

    zplot = pg.zplot
    zpower = pg.zpower
    zenergy = pg.zenergy

    ax.plot(zplot, zenergy, '-')
ax.set_yscale('log')
#%%
temp = pg.get_beamdata('bunching')
ymax = np.max(temp)

fig = plt.figure()
ax1 = fig.add_subplot(5, 1, (1, 3))
ax2 = fig.add_subplot(3, 1, 3)

step = 20
for i in np.arange(0, temp.shape[0], step):
    ax1.plot(temp[i], '-')
    
    ax1.set_title('z = : %.3f m' % zplot[i])
    ax1.set_xlabel('# of slice')
    ax1.set_ylabel('Bunching')
    #ax.set_ylim(0, ymax)
    ax1.grid()
    
    ax2.plot(zplot[:i], ppower[:i]/1e6, '-')
    ax2.set_xlabel(r'$z$ (m)')
    ax2.set_ylabel(r'Total power (MW)')
    ax2.grid()
    
    #display.display(plt.gcf())
    #display.clear_output(wait=True)
    
    plt.pause(0.5)
    plt.cla()
                
    # time.sleep(1)
    # if i<temp.shape[0]-step:
    #     ax1.cla()
    #     ax2.cla()
ax1.set_yscale('log')    
#%% 
x = pg.zplot
current = pg.current

# # version 2
xrms = (pg.get_fielddata('xrms') @ current)/np.sum(current)
yrms = (pg.get_fielddata('yrms') @ current)/np.sum(current)

# version 4

# nn = np.sum(current>0)
# xrms = (pg.get_data('Beam', 'xsize')[:,:nn] @ current[:nn])/np.sum(current)
# yrms = (pg.get_data('Beam', 'ysize')[:,:nn] @ current[:nn])/np.sum(current)
                   
fig, ax = plt.subplots(figsize = (5, 3))

ax.plot(x, xrms*1e3, 'r-')
ax.plot(x, yrms*1e3, 'b-')
ax.grid()

ax.set_xlabel(r'$z$ (m)')
ax.set_ylabel(r'RMS size (mm)')
ax.legend([r'$x$', r'$y$'])

#%% Genesis1.3-version2 modify theta from existing particle output
import h5py

workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim\2022\THzComm2\2nC\SASE\002-G2'
os.chdir(workdir)

fname = 'pithz.out.par'
ALL = np.fromfile(fname, 'float')

Qb = 2e-9
sigma_z = 2e-3
sigma_t = sigma_z/g_c
Ipeak = Qb/np.sqrt(2*np.pi)/sigma_t

lambda0 = 100e-6
k0 = 2*np.pi/lambda0

fG = lambda x, mu = 0, sig = sigma_z:Ipeak*np.exp(-(x-mu)**2/2/sig**2)

c_sig = 2
Ns_b = int(sigma_z*c_sig*2/lambda0)//2*2+2

Ns = 150
Nm_slice = 32768

Nm_slices = 0
stat = []
for islice in np.arange(0, Ns_b+2):
    
    start = islice*Nm_slice*6+Nm_slice
    end   = start+Nm_slice
    theta = ALL[start:end]
    
    x1 = -c_sig*sigma_z+islice*lambda0
    a = (x1-0.5*lambda0)/sigma_z
    b = (x1+0.5*lambda0)/sigma_z
    
    rr = np.array([HaltonNorm2(Nm_slices+j, a = a, b = b) for j in np.arange(Nm_slice)])-a
    
    current1 = fG(x1)
    theta1 = 2*np.pi*(rr*sigma_z/lambda0)-np.pi
    
    Nm_slices += Nm_slice
    
    print(islice, current1)
    
    bf2 = np.abs(np.mean(np.exp(1j*theta)))**2
    bf2_new = np.abs(np.mean(np.exp(1j*theta1)))**2
    stat.append([islice, x1, current1, bf2, bf2_new])

    ALL[start:end] = theta1[:]
    
with open('temp.out.par', 'wb') as f_handle:
    f_handle.write(ALL.tobytes('F'))
    
stat = np.array(stat)

fig, ax = plt.subplots(figsize = (4, 6), nrows = 2, sharex = True)
ax1, ax2 = ax.flatten()
ax1.plot(stat[:,0], stat[:,2])

ax2.plot(stat[:,0], stat[:,3])
ax2.plot(stat[:,0], stat[:,4])

ax2.set_yscale('log')

#%% Visualize the slices from Genesis1.3-version2
# ALL is the data loaded from binary format
fname = 'temp.5.out.par'
ALL = np.fromfile(fname, 'float')

Ns = 150
Nm_slice = 4096*8

fig, ax = plt.subplots(ncols = 3, nrows = 3, figsize = (9, 9))

islice = 40
aslice = ALL[islice*Nm_slice*6:(islice+1)*Nm_slice*6].reshape((6, Nm_slice)).T
gamma, theta, xx, yy, px, py = aslice.T[:]

ax[0,0].hist2d(theta, gamma, bins = 50, cmin = 1e-9)
ax[1,0].hist2d(xx*1e3, px, bins = 50, cmin = 1e-9)
ax[2,0].hist2d(yy*1e3, py, bins = 50, cmin = 1e-9)

for icoor, axis in enumerate(ax[:,1:].flatten()):
    
    cc = 1
    if icoor in [2,3]:
        cc = 1e3
    coor = aslice[:,icoor]*cc
    
    axis.hist(coor, bins = 50, histtype = r'step')
    axis.set_title('RMS = %.2f' % (np.std(coor)), fontsize = 10)

# Bunching along the slices
# Ns = 150
# Nm_slice = 32768

Nm_slices = 0
stat = []
for islice in np.arange(0, Ns_b+2):
    
    start = islice*Nm_slice*6+Nm_slice
    end   = start+Nm_slice
    theta = ALL[start:end]
    
    bf2 = np.abs(np.mean(np.exp(1j*theta)))**2
    
    stat.append([islice, bf2])

stat = np.array(stat)    
with open('bunching@'+fname, 'wb') as f_handle:
    np.savetxt(f_handle, stat, fmt = '%14.6E')

fig, ax1 = plt.subplots(figsize = (5, 4))
ax1.plot(stat[:,0], stat[:,1])
ax1.set_yscale('log')

#%% Genesis1.3-version2 produce slices as input

seed = 5
for seed in [1, 2, 3, 4, 6, 7, 10]:
    beamAstra = '../368A.2809.003'
    #beamGen = '../368A.2809.003.gen'
    
    # x, px, y, py, z, gamma, px=beta_x*gamma
    dist = astra2gen(beamAstra, save = False)
    z = dist[:,4]
    z -= z.mean()
    
    # plt.figure()
    # plt.hist(z*1e3, bins = 50, histtype = r'step')
    
    # Resonant wavelength
    lambda0 = 100e-6
    k0 = 2*np.pi/lambda0
    
    Ns = int((z.max()-z.min())/lambda0)//2*2-4
    Ntail = -Ns//2
    
    # Setup the beam distribution, considering gaussian distribution with same peak current
    Qb = 2e-9
    sigma_z = 2e-3
    sigma_t = sigma_z/g_c
    Ipeak = Qb/np.sqrt(2*np.pi)/sigma_t
    
    znorm = z/sigma_z
    
    # Define the gaussian distribution
    fG = lambda x, mu = 0, sig = sigma_z:Ipeak*np.exp(-(x-mu)**2/2/sig**2)
    
    # Cutoff the gaussian distribution
    c_sig = 2
    c_sig = Ns//2*lambda0/sigma_z
    Ns_b = int(sigma_z*c_sig*2/lambda0)//2*2+2
    
    # Total number of slices and number of macro particles per slice
    Ns = 150
    Nm_slice = 32768
    
    # Initiate the array for storing the coordinates of all slices
    # For Genesis1.3-version2, the coordinates are saved in the order of
    # gamma, theta, x, y, px, py
    
    ALL_slices = np.zeros((Nm_slice*Ns*6))
    
    Nm_slices = 0
    stat = []
    for islice in np.arange(0, Ns_b+2):
        
        # starting and ending indices of one slice
        start = islice*Nm_slice*6
        end   = start+Nm_slice*6
        
        x1 = -c_sig*sigma_z+islice*lambda0
        a = (x1-0.5*lambda0)/sigma_z
        b = (x1+0.5*lambda0)/sigma_z
        
        rr = np.array([HaltonNorm2(Nm_slices+j, a = a, b = b, i = seed) for j in np.arange(Nm_slice)])-a
        
        current1 = fG(x1)
        theta1 = 2*np.pi*(rr*sigma_z/lambda0)-np.pi
        
        Nm_slices += Nm_slice
        
        print(islice, current1)
        
        bf2 = np.abs(np.mean(np.exp(1j*theta)))**2
        bf2_new = np.abs(np.mean(np.exp(1j*theta1)))**2
        stat.append([islice, x1, current1, bf2, bf2_new])
        
        select = (znorm>a-0.02)*(znorm<b+0.02); print(np.sum(select))
        ilast = islice
        if np.sum(select)>0:
            temp, _ = resampleParticles(dist[select], mpart = Nm_slice)
            
            xx = temp[:,0]
            yy = temp[:,2]
            xp = temp[:,1]
            yp = temp[:,3]
            gamma = temp[:,-1]
            
            aslice = np.concatenate((gamma, theta1, xx, yy, xp, yp), axis = 0)
            ALL[start:end] = aslice[:]
        else:
            break
        
    for islice in np.arange(ilast, Ns):
        
        # starting and ending indices of one slice
        start = islice*Nm_slice*6
        end   = start+Nm_slice*6
        
        x1 = -c_sig*sigma_z+islice*lambda0
        a = (x1-0.5*lambda0)/sigma_z
        b = (x1+0.5*lambda0)/sigma_z
        
        rr = np.array([Halton(seed, Nm_slices+j) for j in np.arange(Nm_slice)])
        
        current1 = fG(x1)
        theta1 = 2*np.pi*(rr*sigma_z/lambda0)-np.pi
        
        Nm_slices += Nm_slice
        
        print(islice, current1)
        
        bf2_new = np.abs(np.mean(np.exp(1j*theta1)))**2
        stat.append([islice, x1, current1, bf2, bf2_new])
    
        xx = np.ones(Nm_slice)*dist[:,0].mean()
        yy = np.ones(Nm_slice)*dist[:,2].mean()
        xp = np.ones(Nm_slice)*dist[:,1].mean()
        yp = np.ones(Nm_slice)*dist[:,3].mean()
        gamma = np.ones(Nm_slice)*dist[:,-1].mean()
        
        aslice = np.concatenate((gamma, theta1, xx, yy, xp, yp), axis = 0)
        ALL[start:end] = aslice[:]
        
    with open('temp.%d.out.par' % seed, 'wb') as f_handle:
        f_handle.write(ALL.tobytes('F'))

#%%
def Astra2Genesis2Slices(inputName = None, inputDist = None, outputName = 'temp',
                         seed = 8, npart = 8192*4, nslice = None, lambda0 = 100e-6,
                         Qscale = 1, zscale = 1):
    '''
    Convert `Astra` distribution into slice based distribution for `Genesis1.3 Version2`.

    Parameters
    ----------
    inputName : string, optional
        `Astra` file name. The default is None.
    inputDist : 6D array, optional
        Input distribution from `Astra`. The default is None. One of `inputName` 
        or `inputDist` must be defined.
    outputName : string, optional
        File name of the Genesis distribution. The default is `temp`.
    seed : int, optional
        Seeding number of the Hammersley sequence. The default is 8.
    npart : int, optional
        Number of macro particles of each slice. The default is 8192*4.
    nslice : int, optional
        Number of slices. The default is None.
    lambda0 : double, optional
        Resonanat wavelength. The default is 100e-6.
    Qscale : double, optional
        Scale factor of the bunch charge from `Astra` distribution. The default is 1.
    zscale : double, optional
        Scale factor of the bunch length from `Astra` distribution. The default is 1.

    Returns
    -------
    Ntail : int
        Same as `Ntail` in `Genesis1.3`.
    curpeak : float
        Peak current, A.
    curlen : float
        RMS length of the current profile, positive for Gaussian and negative for flattop,
        meter.

    '''
    # dist: x, px, y, py, z, gamma. px=beta_x*gamma
    #dist = astra2gen(beamAstra, save = False)
    
    if inputName != None:   
        dist = pd_loadtxt(inputName)
        dist[1:,2] += dist[0,2]
        dist[1:,5] += dist[0,5]
    elif inputDist != None:
        dist = inputDist[:,0:6]
    else:
        print('No input file or data!')
        return
        
    if outputName == None:
        outputName = 'temp'
    outputName += str.format('.%d.out.par' % seed)
    print('The distribution is saved to '+outputName)
    
    x, y, z = dist[:,0:3].T[:]
    z -= z.min()
    
    zmin, zmax = z.min()*1e3, z.max()*1e3
    print('Bunch length in mm: ', zmax-zmin)
    
    Px, Py, Pz = dist[:,3:6].T[:]
    
    px = Px/g_mec2/1e6
    py = Py/g_mec2/1e6
    gg = np.sqrt(1+(Px**2+Py**2+Pz**2)/g_mec2/1e6/g_mec2/1e6) # gamma
    
    newdist = np.zeros(dist[:,0:6].shape)
    newdist[:,0] = x
    newdist[:,1] = px[:]
    newdist[:,2] = y
    newdist[:,3] = py
    newdist[:,4] = z
    newdist[:,5] = gg
    
    Qb = np.sum(dist[:,-3])*1e-9*Qscale
    Qb = np.abs(Qb)
    
    z = newdist[:,4]
    z -= z.mean()
    z *= zscale
    
    # Resonant wavelength
    #lambda0 = 100e-6
    #lambda0 = xlamb
    k0 = 2*np.pi/lambda0
    
    # Calculate Ntail = number of slices of the e bunch/2
    Ns_bunch = int((z.max()-z.min())/lambda0)//2*2-4
    Ntail = -Ns_bunch//2
    
    # Setup the beam distribution, considering gaussian distribution with same peak current
    #Qb = 2e-9
    
    #sigma_z = 2e-3
    sigma_z = z.std()
    
    sigma_t = sigma_z/g_c
    Ipeak = Qb/np.sqrt(2*np.pi)/sigma_t
    
    curpeak, curlen = Ipeak, sigma_z
    
    znorm = z/sigma_z
    
    # Define the gaussian distribution
    fG = lambda x, mu = 0, sig = sigma_z:Ipeak*np.exp(-(x-mu)**2/2/sig**2)
    
    # Cutoff the gaussian distribution
    c_sig = 3
    c_sig = Ns_bunch//2*lambda0/sigma_z
    Ns_b = int(sigma_z*c_sig*2/lambda0)//2*2+2
    
    # Total number of slices and number of macro particles per slice
    #Ns = 150
    if nslice == None:
        Ns = Ns_bunch+120
    else:
        Ns = nslice
        
    #Nm_slice = 32768
    Nm_slice = npart
    
    # Initiate the array for storing the coordinates of all slices
    # For Genesis1.3-version2, the coordinates are saved in the order of
    # gamma, theta, x, y, px, py
    
    ALL = np.zeros((Nm_slice*Ns*6))
    
    rn = HaltonN(Nm_slice*(Ns_b+2), iprime = seed)
    
    Nm_slices = 0
    stat = []
    for islice in np.arange(0, Ns_b+2):
        
        # starting and ending indices of one slice
        start = islice*Nm_slice*6
        end   = start+Nm_slice*6
        
        x1 = -c_sig*sigma_z+islice*lambda0
        a = (x1-0.5*lambda0)/sigma_z
        b = (x1+0.5*lambda0)/sigma_z
        
        #rr = np.array([HaltonNorm2(Nm_slices+j, a = a, b = b, i = seed) 
        #               for j in np.arange(Nm_slice)])-a
        
        samples = rn[Nm_slices:Nm_slices+Nm_slice]
        rr = HaltonNNorm2(samples, a = a, b = b)-a
        
        current1 = fG(x1)
        theta1 = 2*np.pi*(rr*sigma_z/lambda0)-np.pi
        
        Nm_slices += Nm_slice
        
        print('# of slices: ', islice, ', current = ', current1, ' A')
        
        #bf2 = np.abs(np.mean(np.exp(1j*theta)))**2
        bf2 = np.abs(np.mean(np.exp(1j*theta1)))**2
        stat.append([islice, x1, current1, bf2])
        
        select = (znorm>a-0.02)*(znorm<b+0.02); #print(np.sum(select))
        ilast = islice
        if np.sum(select)>0:
            temp, _ = resampleParticles(newdist[select], mpart = Nm_slice)
            
            xx = temp[:,0]
            yy = temp[:,2]
            xp = temp[:,1]
            yp = temp[:,3]
            gamma = temp[:,-1]
            
            aslice = np.concatenate((gamma, theta1, xx, yy, xp, yp), axis = 0)
            ALL[start:end] = aslice[:]
        else:
            break
        
    for islice in np.arange(ilast, Ns):
        
        # starting and ending indices of one slice
        start = islice*Nm_slice*6
        end   = start+Nm_slice*6
        
        x1 = -c_sig*sigma_z+islice*lambda0
        a = (x1-0.5*lambda0)/sigma_z
        b = (x1+0.5*lambda0)/sigma_z
        
        #rr = np.array([Halton(seed, Nm_slices+j) for j in np.arange(Nm_slice)])
        rr = HaltonN(Nm_slice, iprime = seed)
        
        current1 = fG(x1)*0
        theta1 = 2*np.pi*(rr*sigma_z/lambda0)-np.pi
        
        Nm_slices += Nm_slice
        
        print('# of slices: ', islice, ', current = ', current1, ' A')
        
        bf2 = np.abs(np.mean(np.exp(1j*theta1)))**2
        stat.append([islice, x1, current1, bf2])
    
        xx = np.ones(Nm_slice)*newdist[:,0].mean()
        yy = np.ones(Nm_slice)*newdist[:,2].mean()
        xp = np.ones(Nm_slice)*newdist[:,1].mean()
        yp = np.ones(Nm_slice)*newdist[:,3].mean()
        gamma = np.ones(Nm_slice)*newdist[:,-1].mean()
        
        aslice = np.concatenate((gamma, theta1, xx, yy, xp, yp), axis = 0)
        ALL[start:end] = aslice[:]
        
    with open(outputName, 'wb') as f_handle:
        f_handle.write(ALL.tobytes('F'))
    
    stat = np.array(stat)
    header = str.format('%12s%14s%14s%14s' % ('islice', 'z', 'I', 'bf2'))
    with open('bunching@'+outputName, 'wb') as f_handle:
        np.savetxt(f_handle, stat, fmt = '%14.6E', 
                   header = header)
    
    fig, ax = plt.subplots(nrows = 2, figsize = (5, 6))
    ax1, ax2 = ax.flatten()
    ax1.plot(stat[:,0], stat[:,-1], '-')
    ax1.set_yscale('log')
    ax1.grid()
    ax1.set_xlabel(r'Number of slices')
    ax1.set_ylabel(r'$|b_f|^2$')
    
    ax2.plot(stat[:,1]*1e3, stat[:,2], '-')
    ax2.grid()
    ax2.set_xlabel(r'$z$ (mm)')
    ax2.set_ylabel(r'$I_{peak}$ (A)')
    
    fig.tight_layout()
    fig.savefig('bunching@'+outputName+'.png')
    
    return curpeak, curlen, Ntail

#%%
workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim3\2023\G2-demo'
os.chdir(workdir)

fname = '368A.2809.003'
curpeak, curlen, Ntail = Astra2Genesis2Slices(fname, outputName = 'test', seed = 1, npart = 4096*8)
curpeak, curlen, Ntail = Astra2Genesis2Slices(fname, outputName = 'test', seed = 2, npart = 4096*8)
curpeak, curlen, Ntail = Astra2Genesis2Slices(fname, outputName = 'test', seed = 8, npart = 4096*8)
curpeak, curlen, Ntail = Astra2Genesis2Slices(fname, outputName = 'test', seed = 9, 
                                              npart = 4096*8)

#%%
outputName = 'temp.8.out.par'
stat = np.loadtxt('bunching@'+outputName)

fig, ax = plt.subplots(nrows = 2, figsize = (5, 6))
ax1, ax2 = ax.flatten()
ax1.plot(stat[:,0], stat[:,-1], '-')
ax1.set_yscale('log')
ax1.grid()
ax1.set_xlabel(r'Number of slices')
ax1.set_ylabel(r'$|b_f|^2$')

ax2.plot(stat[:,1]*1e3, stat[:,2], '-')
ax2.grid()
ax2.set_xlabel(r'$z$ (mm)')
ax2.set_ylabel(r'$I_{peak}$ (A)')


fig.tight_layout()

fig.savefig('bunching@'+outputName+'.png')


#%%
workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim3\2023\THUtest'
os.chdir(workdir)

dist = np.loadtxt('TTX.beam', skiprows = 4)

fig, ax = plt.subplots()
ax.hist2d(dist[:,4], dist[:,5], bins = 100, cmin = 1e-9)

#%%
#fig, ax = plt.subplots()
#ax.hist(dist[:,4]/g_c*1e12, bins = 100, histtype = r'step')

fname = 'TTX.laser'
dist = np.loadtxt('TTX.laser', skiprows = 3)

fig, ax = plt.subplots()
ax.plot(dist[:,0], dist[:,1]*100e-6/g_c)
ax.set_yscale('log')

