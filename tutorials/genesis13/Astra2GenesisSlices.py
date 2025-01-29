import sys
sys.path.append('/afs/ifh.de/group/pitz/data/lixiangk/work/sync/python')
#sys.path.append('\\afs\ifh.de\group\pitz\data\lixiangk\work\sync\python')

from timeit import default_timer
import numpy as np

import time

# def CDFFlattop(z, Lz = 1, rz = 0.1):
#     r = -rz*z
#% Define the Hammersley sequence and the function to convert Astra distrution 
#   to slices for Genesis1.3 Version 2
from xkinterface.interface import *

from scipy.stats import qmc
erf = sp.special.erf
erfinv = sp.special.erfinv

#%% Define the functions for smoothing and calculating the density of the current profile
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

def HaltonNNorm2(Hsamples, a = -10, b = 10):
    '''
    Generate samples following normal distribution with random
    samples generated from Hammersley sequence.
    sequence 

    Parameters
    ----------
    Hsamples : 1D
        Random samples (0, 1) generated from Hammersley sequence.
    a : float, optional
        Left position. The default is -10.
    b : float, optional
        Right position. The default is 10.
        
    Returns
    -------
    samples : 1D
        Samples following normal distribution.
        
    '''
    # F(x, a) = (erf(x/np.sqrt(2))-erf(a/np.sqrt(2)))/2/ccc -> 0, 1
    Fa = erf(a/np.sqrt(2))
    Fb = erf(b/np.sqrt(2))
    
    #r1 = HaltonN(n, iprime)
    #r1 = Hsamples
    samples = erfinv(2*Hsamples*(Fb-Fa)/2+Fa)*np.sqrt(2)
    
    return samples

def InvCDF(samples, wbin = 100e-6, degree = 4):
    '''
    Calculate the inverse CDF of given samples

    Parameters
    ----------
    samples : 1D
        Samples for calculating the inverse CDF.

    Returns
    -------
    Fx_inv : function
        The inverse CDF.
    '''
    smin = samples.min()
    smax = samples.max()
    bins = int((smax-smin)/wbin)*1;  #print(bins)
    
    weights, edges = np.histogram(samples, bins = bins)
    ### Add zeros to avoid abrupt change of density near tails
    weights = np.concatenate(([0]*9, weights, [0]*9))
    ds = edges[1]-edges[0]
    
    # edges = np.concatenate((np.array([-3, -2, -1])*ds+edges[0], 
    #                         edges, 
    #                         np.array([1, 2, 3])*ds+edges[-1]))
    ###
    # ### Smooth the current profile, has to be odd number
    # def easy_smooth(y, box_pts):
    #     '''
    #     Smooth the data points using a moving box
    #     Parameters
    #       y: 1-D array to be smoothed
    #       box_pts: number of points in the moving box
    #     Returns
    #       y_smooth: smoothed 1-D array
    #     '''
    #     box = np.ones(box_pts)/box_pts
    #     y_smooth = np.convolve(y, box, mode='same')
    #     return y_smooth
    # #weights = easy_smooth(weights, 32+1)
    
    from scipy.signal import savgol_filter
    
    window_size = 16+1
    #degree = 4  # Degree of the polynomial to fit
    plt.figure(); plt.plot(edges[:], weights[8:-9],label='ini')
    
    weights = savgol_filter(weights, window_size, degree)
    
    weights[weights<10] = 0; #weights = weights[3:-3]
    plt.plot(edges[:], weights[8:-9], '-',label='fit')
    plt.savefig('Current_profile_smoothed.png')
    ###
    
    edges, weights = edges[:]+ds/2, weights[8:-9]
    ss = weights>0
    edges = edges[ss]
    weights = weights[ss]
    
    
    weights[0] = 0
    weights = weights/np.sum(weights)
    
    fx = interp1d(edges, weights, kind = 'linear', 
                  bounds_error = False, fill_value = 0)
    
    for j in np.arange(1, len(weights)):
        weights[j] += weights[j-1]
    #weights = np.concatenate(([0], weights))
    #weights = weights/weights[-1]
    ###
    
    Fx = interp1d(edges, weights, kind = 'cubic', 
                  bounds_error = False, fill_value = (0, 1))
    Fx_inv = interp1d(weights, edges, kind = 'cubic',
                      bounds_error = False, fill_value = 0)
    
    return Fx_inv, Fx, fx

def HaltonNUser(Hsamples, a = -10, b = 10, Fx_inv = None, Fx = None):
    '''
    Generate samples following normal distribution with random
    samples generated from Hammersley sequence.
    sequence 
    
    Parameters
    ----------
    Hsamples : 1D
        Random samples (0, 1) generated from Hammersley sequence.
    a : float, optional
        Left position. The default is -10.
    b : float, optional
        Right position. The default is 10.

    Returns
    -------
    samples : 1D
        Samples following normal distribution.
        
    '''
    # F(x, a) = (erf(x/np.sqrt(2))-erf(a/np.sqrt(2)))/2/ccc -> 0, 1
    Fa = Fx(a)
    Fb = Fx(b)
    
    Hsamples = Hsamples*(Fb-Fa)+Fa
    samples = Fx_inv(Hsamples)
    
    return samples

# #%%
# import os

# workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim3\2023\THzComm2'
# os.chdir(workdir)

# fname = '368A.2809.002.1000'
# data = pd_loadtxt(fname)

# samples = np.copy(data[1:,2])-np.mean(data[1:,2])
# samples *= 1e3 # mm

# #%%
# fig, ax = plt.subplots(figsize = (5, 5))
# ax.hist(samples, bins = 100, histtype = r'step', range = (-5, 5),
#         cumulative = False)#, density = True)

# Fx_inv, Fx, fx = InvCDF(samples*1e-3)

# xx = np.linspace(-5, 5, 1001)
# yy = fx(xx*1e-3)

# ax.plot(xx, yy, '-*')

# xx = np.linspace(-5, 5, 101)
# yy = fx(xx*1e-3)

# ax.plot(xx, yy, '-')

#%% Astra to Genesis1.3 Version 2
def Astra2Genesis2Slices(inputName = None, inputDist = None, outputName = 'temp',
                         seed = 8, npart = 8192*4, nslice = None, lambda0 = 100e-6,
                         Qscale = 1, zscale = 1, bunch = 0):
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
    Ns: int
        Same as `nslice` in Genesis1.3 V2
    Ntail: int
        Same as `ntail` in Genesis1.3 V2
    outputName: string
        Name of the file where all slices are saved, to be used as `partfile` in Genesis1.3 V2

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
    
    # Resonant wavelength
    #lambda0 = 100e-6
    k0 = 2*np.pi/lambda0
    
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
    
    # RMS length after scaling
    sigma_z = z.std()
    sigma_t = sigma_z/g_c
    
    Ipeak = Qb/np.sqrt(2*np.pi)/sigma_t
    znorm = z/sigma_z
    
    curpeak, curlen = Ipeak, sigma_z
    
    # Define the gaussian distribution
    fG = lambda x, mu = 0, sig = sigma_z:Ipeak*np.exp(-(x-mu)**2/2/sig**2)
    
    # Calculate Ntail = number of slices of the e bunch/2, based on the full length
    Ns_b = int((z.max()-z.min())/lambda0)//2*2-4
    
    # Or Cutoff the gaussian distribution
    c_sig = 2.5
    c_sig = Ns_b//2*lambda0/sigma_z
    Ns_b = int(sigma_z*c_sig*2/lambda0)//2*2+2
    
    Ntail = -Ns_b//2
    
    # Total number of slices and number of macro particles per slice
    #Ns = 150
    if nslice == None:
        Ns = Ns_b+100
    else:
        Ns = nslice
        
    #Nm_slice = 32768
    Nm_slice = npart
    
    # Initiate the array for storing the coordinates of all slices
    # For Genesis1.3-version2, the coordinates are saved in the order of
    # gamma, theta, x, y, px, py
    ALL = np.zeros((Nm_slice*Ns*6))
    
    # Generate Hammersley sequence
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
        #samples = np.random.rand(Nm_slice)
        rr = HaltonNNorm2(samples, a = a, b = b)-a
        
        current1 = fG(x1)
        theta1 = 2*np.pi*(rr*sigma_z/lambda0)-np.pi
        if bunch>0: # 𝜃:=𝜃−2𝑏 sin⁡〖(𝜃−𝜙)〗
            theta1 -= 2*bunch*np.sin(theta1)
            
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
    
    select = (znorm>-0.1)*(znorm<0.1)
    temp, _ = resampleParticles(newdist[select], mpart = Nm_slice)
    
    for islice in np.arange(ilast, Ns):
        
        # starting and ending indices of one slice
        start = islice*Nm_slice*6
        end   = start+Nm_slice*6
        
        x1 = -c_sig*sigma_z+islice*lambda0
        
        #rr = np.array([Halton(seed, Nm_slices+j) for j in np.arange(Nm_slice)])
        rr = HaltonN(Nm_slice, iprime = seed)
        
        current1 = fG(x1)*0
        theta1 = 2*np.pi*(rr*sigma_z/lambda0)-np.pi
        
        Nm_slices += Nm_slice
        
        print('# of slices: ', islice, ', current = ', current1, ' A')
        
        bf2 = np.abs(np.mean(np.exp(1j*theta1)))**2
        stat.append([islice, x1, current1, bf2])
        
        xx = temp[:,0]
        yy = temp[:,2]
        xp = temp[:,1]
        yp = temp[:,3]
        gamma = temp[:,-1]
            
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
    
    return curpeak, curlen, Ns, Ntail, outputName

#%% Astra to Warp with Hamsley sequence in z
def Astra2WarpHaltonZ(inputName = None, inputDist = None, 
                         outputName = 'temp',
                         seed = 8, lambda0 = 100e-6,
                         Qscale = 1, zscale = 1, bunch = 0,
                         useHammersley = 1):
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
    lambda0 : double, optional
        Resonanat wavelength. The default is 100e-6.
    Qscale : double, optional
        Scale factor of the bunch charge from `Astra` distribution. The default is 1.
    zscale : double, optional
        Scale factor of the bunch length from `Astra` distribution. The default is 1.

    Returns
    -------
    outputName: string
        Name of the file where all slices are saved, to be used as `partfile` in Genesis1.3 V2

    '''
    # dist: x, px, y, py, z, gamma. px=beta_x*gamma
    #dist = astra2gen(beamAstra, save = False)
    astra2warp
    if inputName != None:   
        dist = pd_loadtxt(inputName)
        dist[1:,2] += dist[0,2]
        dist[1:,5] += dist[0,5]
        #select = dist[:,-1]>0
        #dist = dist[select]
        
    elif inputDist != None:
        dist = inputDist
        
    else:
        print('No input file or data!')
        return
        
    if outputName == None:
        outputName = 'temp'
    outputName += str.format('.%d.warp' % seed)
    print('The distribution is saved to '+outputName)
    
    # Resonant wavelength
    #lambda0 = 100e-6
    k0 = 2*np.pi/lambda0
    
    x, y, z = dist[:,0:3].T[:]
    z -= z.min()
    
    zmin, zmax = z.min()*1e3, z.max()*1e3
    print('Bunch length in mm: ', zmax-zmin)
    
    select = (dist[:,9]>0); dist = dist[select]
    dist[:,3] = dist[:,3]/g_mec2/1e6 # beta_x*gamma
    dist[:,4] = dist[:,4]/g_mec2/1e6 # beta_y*gamma
    dist[:,5] = dist[:,5]/g_mec2/1e6 # beta_z*gamma
    
    Q_coef = -1
    dist[:,6] = dist[:,7]*1e-9/g_qe*Q_coef # number of electrons for each macro particle
    Qb = np.sum(dist[:,6])*g_qe
    
    z = dist[:,2]
    z -= z.mean()
    z *= zscale
    
    zmin = z.min()
    zmax = z.max()
    
    # RMS length after scaling
    sigma_z = z.std()
    sigma_t = sigma_z/g_c
    
    # Generate Hammersley sequence
    Nm = len(dist)
    Ne = int(Qb/g_qe)
    
    rn = HaltonN(Nm, iprime = seed)
    
    wbin = lambda0*zscale/1
    #nbin = lambda0/wbin
    #nbin = 4
    Fx_inv, Fx, fx = InvCDF(z, wbin = wbin)
    
    a = int(zmin/lambda0)*lambda0
    b = int(zmax/lambda0)*lambda0
    
    #samples = rn[Nm_slices:Nm_slices+Nm_slice]
    #rr = HaltonNNorm2(samples, a = a, b = b)-a
    
    # temp, _ = resampleParticles(newdist[select], 
    #                             mpart = npart)
    
    eps_scale = 1 #np.sqrt(2)
    dist[:,0] *= eps_scale
    dist[:,1] *= eps_scale
    dist[:,3] *= eps_scale
    dist[:,4] *= eps_scale
    
    z1 = dist[:,2]
    #gamma = dist[:,-1]
    
    ### Set all gamma the same, just for test
    # g0 = np.mean(gamma)
    # gamma = g0*np.ones(len(gamma))
    ###
    
    # Generate samples between 0 and 1
    if useHammersley: 
        samples = rn
    else:
        r0 = (np.arange(1, Nm+1)-0.5)/Nm
        delta = np.sqrt(3*Nm/Ne)
        
        rr = np.random.rand(Nm)
        dr = ((2*rr-1)*delta/2/np.pi)*1
        
        samples = r0+dr
        
    rr = HaltonNUser(samples, a = a, b = b, 
                     Fx_inv = Fx_inv, Fx = Fx)
    
    ### this will disable the current file, just for test
    # rr = a+(b-a)*samples # 
    
    # theta1_ = 2*np.pi*((rr-a)/lambda0)-np.pi
    z1_ = rr
    
    arg = np.argsort(z1)
    z1[arg] = np.sort(z1_)
    
    #np.savetxt(outputName, dist[:,0:7], fmt = '%20.12E')
    
    gamma = np.sqrt(1+dist[:,3]**2+dist[:,4]**2+dist[:,5]**2)
    vz = dist[:,5]/gamma*g_c
    
    vzmean = np.mean(vz)
    
    zz = np.copy(z1)+(vz-vzmean)*1000e-12
    
    islice = 0
    z_ = a
    stat = []
    while z_ <= b:
        zmin_ = z_-lambda0/2
        zmax_ = z_+lambda0/2
        
        select_ = (zz>zmin_)*(zz<zmax_)
        
        zz_ = zz[select_]
        
        theta = (zz_-zmin_)/lambda0*2*np.pi
        bf = np.mean(np.exp(1j*theta))
        
        current1 = (Fx(zmax_)-Fx(zmin_))*Qb/(lambda0/g_c)
        Ne_slice = (Fx(zmax_)-Fx(zmin_))*Ne
        
        bf2 = np.abs(bf)**2
        phase = np.angle(bf)
        stat.append([islice, z_, current1, bf2, bf2*(current1)**2, Ne_slice, 
                     phase, np.real(bf), np.imag(bf)])
    
        islice += 1
        z_ += lambda0
    
    stat = np.array(stat); curpeak = np.max(stat[:,2])
    header = str.format('%12s%14s%14s%14s' % ('islice', 'z', 'I', 'bf2'))
    with open('bunching_%d@' % seed + outputName, 'wb') as f_handle:
        np.savetxt(f_handle, stat, fmt = '%14.6E', 
                   header = header)
    
    
    ilast = -1
    #fig, ax = plt.subplots(nrows = 4, figsize = (5, 9), sharex = True)
    ax1, ax2, ax3, ax4 = ax.flatten()
    ax1.plot(stat[:ilast,0], stat[:ilast,3], '-')
    ax1.plot(stat[:ilast,0], 1/stat[:ilast,5], 'k--')
    ax1.set_yscale('log')
    ax1.grid()
    #ax1.set_xlabel(r'Number of slices')
    ax1.set_ylabel(r'$|b_f|^2$')
    ax1.set_ylim(1e-10, 1)
    
    # ax2.plot(stat[:ilast,0], stat[:ilast,-3], '-'); 
    # print('Sum of bunching: ', np.sum(stat[:ilast,-3]))
    # ax2.set_yscale('log')
    # ax2.set_ylabel(r'$N^2|b_f|^2$')
    
    #ax2.plot(stat[:ilast,0], stat[:ilast,8]/stat[:ilast,7], '-')
    ax2.plot(stat[:ilast,0], stat[:ilast,6], '-')
    ax2.set_ylabel('Phase (rad)')
    
    ax3.plot(stat[:ilast,0], stat[:ilast,7], '-')
    ax3.plot(stat[:ilast,0], stat[:ilast,8], '-')
    
    #ax4.plot(stat[:ilast,1]*1e3, stat[:ilast,2], '-')
    ax4.plot(stat[:ilast,0], stat[:ilast,2], '-')
    ax4.grid()
    ax4.set_xlabel(r'$z$ (mm)')
    ax4.set_ylabel(r'$I_{peak}$ (A)')
    ax4.set_xlabel(r'Number of slices')
    
    fig.tight_layout()
    fig.savefig('bunching_%d@' % seed + outputName+'.png')
    
    return outputName

#%% Making bunching plot from slice-wise statistics
def bunchingPlot(stat, ilast = -1):
    
    
    fig, ax = plt.subplots(nrows = 4, figsize = (5, 9), sharex = True)
    ax1, ax2, ax3, ax4 = ax.flatten()
    ax1.plot(stat[:ilast,0], stat[:ilast,3], '-')
    ax1.plot(stat[:ilast,0], 1/stat[:ilast,5], 'k--')
    ax1.set_yscale('log')
    ax1.grid()
    #ax1.set_xlabel(r'Number of slices')
    ax1.set_ylabel(r'$|b_f|^2$')
    ax1.set_ylim(1e-10, 1)
    
    # ax2.plot(stat[:ilast,0], stat[:ilast,-3], '-'); 
    # print('Sum of bunching: ', np.sum(stat[:ilast,-3]))
    # ax2.set_yscale('log')
    # ax2.set_ylabel(r'$N^2|b_f|^2$')
    
    #ax2.plot(stat[:ilast,0], stat[:ilast,8]/stat[:ilast,7], '-')
    ax2.plot(stat[:ilast,0], stat[:ilast,6], '-')
    ax2.set_ylabel('Phase (rad)')
    
    ax3.plot(stat[:ilast,0], stat[:ilast,7], '-')
    ax3.plot(stat[:ilast,0], stat[:ilast,8], '-')
    
    #ax4.plot(stat[:ilast,1]*1e3, stat[:ilast,2], '-')
    ax4.plot(stat[:ilast,0], stat[:ilast,2], '-')
    ax4.grid()
    ax4.set_xlabel(r'$z$ (mm)')
    ax4.set_ylabel(r'$I_{peak}$ (A)')
    ax4.set_xlabel(r'Number of slices')
    
    fig.tight_layout()
    fig.savefig('bunching@'+outputName+'.png')

#%% Astra to Genesis1.3 Version 3 or 4
def Astra2GenesisSlices(inputName = None, inputDist = None, 
                         outputName = 'temp',
                         seed = 8, npart = 8192*4, nslice = None, 
                         lambda0 = 100e-6, degree = 4, nperlambda = 1,
                         Qscale = 1, zscale = 1, bunch = 0,
                         useHammersley = 1,
                         version = 4, **kwargs):
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
    bunch : float, optional
        Bunching factor introduced evenly to each slice. The default is zero.
    useHammersley: int, optional
        If 1, random samples will be generated using Hammersley sequence, otherwise 
        using the method proposed by C. Penman, 1992. The default is 1.
    version: int, optional
        Genesis1.3 version. The default is 4.
        
    Returns
    -------
    Ntail : int
        Same as `Ntail` in `Genesis1.3`.
    curpeak : float
        Peak current, A.
    curlen : float
        RMS length of the current profile, positive for Gaussian and negative for flattop,
        meter.
    Ns: int
        Same as `nslice` in Genesis1.3 V2
    Ntail: int
        Same as `ntail` in Genesis1.3 V2
    outputName: string
        Name of the file where all slices are saved, to be used as `partfile` in Genesis1.3 V2

    '''
    # dist: x, px, y, py, z, gamma. px=beta_x*gamma
    #dist = astra2gen(beamAstra, save = False)
    
    if inputName != None:   
        dist = pd_loadtxt(inputName)
        dist[1:,2] += dist[0,2]
        dist[1:,5] += dist[0,5]
        
        beam_P0 = dist[0,5]/1e6
    elif inputDist != None:
        dist = inputDist[:,0:6]
    else:
        print('No input file or data!')
        return
        
    if outputName == None:
        outputName = 'temp'
    outputName += str.format('.%d.out.par.h5' % seed)
    # print('The distribution is saved to '+outputName)
    
    # Resonant wavelength
    #lambda0 = 100e-6
    k0 = 2*np.pi/lambda0
    
    x, y, z = dist[:,0:3].T[:]
    z -= z.min()
    
    zmin, zmax = z.min()*1e3, z.max()*1e3
    # print('Bunch length in mm: ', zmax-zmin)
    
    Px, Py, Pz = dist[:,3:6].T[:]
    
    px = Px/g_mec2/1e6
    py = Py/g_mec2/1e6
    gg = np.sqrt(1+(Px**2+Py**2+Pz**2)/g_mec2/1e6/g_mec2/1e6) # gamma
    
    newdist = np.zeros(dist[:,0:6].shape)
    newdist[:,0] = x
    newdist[:,1] = px[:] # beta*gamma
    newdist[:,2] = y
    newdist[:,3] = py
    newdist[:,4] = z
    newdist[:,5] = gg
    
    Qb = np.sum(dist[:,-3])*1e-9*Qscale
    Qb = np.abs(Qb)
    
    z = newdist[:,4]
    z -= z.mean()
    z *= zscale
    
    #select = (z>-0.75e-3)*(z<0.75e-3)
    #z = z[select]
    #newdist = newdist[select]
    
    zmin = z.min()
    zmax = z.max()
    
    # RMS length after scaling
    sigma_z = z.std()
    sigma_t = sigma_z/g_c
    
    Ipeak = Qb/np.sqrt(2*np.pi)/sigma_t
    znorm = z/sigma_z
    
    curpeak, curlen = Ipeak, sigma_z
    
    # Define the gaussian distribution
    fG = lambda x, mu = 0, sig = sigma_z:Ipeak*np.exp(-(x-mu)**2/2/sig**2)
    
    # Calculate Ntail = number of slices of the e bunch/2, based on the full length
    Ns_b = int((z.max()-z.min())/lambda0)//2*2+1
    
    # Or Cutoff the gaussian distribution
    c_sig = Ns_b//2*lambda0/sigma_z
    
    #c_sig = 2.5
    #Ns_b = int(sigma_z*c_sig*2/lambda0)//2*2+2
    
    Ntail = -Ns_b//2
    
    # Total number of slices and number of macro particles per slice
    if nslice == None:
        nproc = 16
        # Ns_b = (Ns_b//nproc+1)*nproc 
        Ns = Ns_b +100
        Ns = (Ns//nproc+1)*nproc 
        
    else:
        Ns = nslice
        
    Nm_slice = npart
    
    # Generate Hammersley sequence
    rn = HaltonN(Nm_slice*(Ns_b+2), iprime = seed)
    np.random.seed(seed = seed*1234567)
    
    wbin = 100e-6*zscale/nperlambda
    Fx_inv, Fx, fx = InvCDF(z, wbin = wbin, degree = degree)
    
    # saving to hdf5
    f = h5py.File(outputName, 'w')
    
    Nm_slices = 0
    stat = []
    for islice in np.arange(0, Ns_b+1):
        
        # Center of the slice
        x1 = zmin+0*lambda0+islice*lambda0
        
        a = (x1-0.5*lambda0)
        b = (x1+0.5*lambda0)
        
        w1 = 1 # Select particles which are within +/- w1*lambda/2
        select = (z>x1-w1*lambda0/2)*(z<=x1+w1*lambda0/2)
        current0 = np.sum(select)/len(z)*Qb/(lambda0/g_c)/w1
        current1 = (Fx(b)-Fx(a))*Qb/(lambda0/g_c)
        
        Ne_slice = current1*lambda0/g_c/g_qe
        
        # print('# of slices: ', islice, ', current = ', current0, ' / ', current1, ' A')
        
        # For dummy particles (current1 = 0)
        ccc0 = 1
        temp0, _ = resampleParticles(newdist[::11], mpart = Nm_slice)
        
        #select = (znorm>a-0.02)*(znorm<b+0.02); #print(np.sum(select))
        ilast = islice
        if np.sum(select)>0 and current1>0:
            temp, _ = resampleParticles(newdist[select], 
                                        mpart = Nm_slice)
            Nm_slices += Nm_slice
            
            eps_scale = 1 #np.sqrt(2)
            xx = temp[:,0]*eps_scale
            yy = temp[:,2]*eps_scale
            xp = temp[:,1]*eps_scale
            yp = temp[:,3]*eps_scale
            
            theta1 = (temp[:,4]-x1)/lambda0*2*np.pi
            gamma = temp[:,-1]
            
            ### Set all gamma the same, just for test
            # g0 = np.mean(gamma)
            # gamma = g0*np.ones(len(gamma))
            ###
            
            # Generate samples between 0 and 1
            if useHammersley: 
                samples = rn[Nm_slices:Nm_slices+Nm_slice]
            else:
                r0 = (np.arange(1, Nm_slice+1)-0.5)/Nm_slice
                delta = np.sqrt(3*Nm_slice/Ne_slice)
                
                rr = np.random.rand(Nm_slice)
                dr = ((2*rr-1)*delta/2/np.pi)*1
                
                samples = r0+dr
                
            rr = HaltonNUser(samples, a = a, b = b, 
                             Fx_inv = Fx_inv, Fx = Fx)
            
            ### this will disable the current file, just for test
            # rr = a+(b-a)*samples # 
            
            theta1_ = 2*np.pi*((rr-a)/lambda0)-np.pi
                        
            arg = np.argsort(theta1)
            theta1[arg] = np.sort(theta1_)
            ###theta1 = np.copy(theta1_)
            
            if bunch>0: # 𝜃:=𝜃−2𝑏 sin⁡〖(𝜃−𝜙)〗
                theta1 -= 2*bunch*np.sin(theta1+np.pi/2)
            
        else:
            
            xx = temp0[:,0]*ccc0
            yy = temp0[:,2]*ccc0
            xp = temp0[:,1]*ccc0
            yp = temp0[:,3]*ccc0
            #theta1 = np.linspace(-np.pi, np.pi, Nm_slice)
            gamma = temp0[:,-1]
            
            rr = HaltonN(Nm_slice, iprime = seed)
            theta1 = (2*rr-1)*np.pi
            
            current1 = 0
            
            Nm_slices += Nm_slice
            
        slicename = "/slice%06d/" % (islice+1)
        f.create_dataset(slicename+"current", data = [current1])
        f.create_dataset(slicename+"gamma", data = gamma)
        f.create_dataset(slicename+"theta", data = theta1)
        f.create_dataset(slicename+"x", data = xx)
        f.create_dataset(slicename+"y", data = yy)
        f.create_dataset(slicename+"px", data = xp)
        f.create_dataset(slicename+"py", data = yp)
        
        bf = np.mean(np.exp(1j*theta1))
        bf2 = np.abs(bf)**2
        phase = np.angle(bf)
        stat.append([islice, x1, current1, bf2, bf2*(current1)**2, Ne_slice, 
                     phase, np.real(bf), np.imag(bf)])
            
    #select = (znorm>-0.1)*(znorm<0.1)
    for islice in np.arange(ilast+1, Ns):
        
        # starting and ending indices of one slice
        start = islice*Nm_slice*6
        end   = start+Nm_slice*6
        
        x1 = -c_sig*sigma_z+islice*lambda0
        x1 = zmin+0*lambda0+islice*lambda0
        
        #rr = np.array([Halton(seed, Nm_slices+j) for j in np.arange(Nm_slice)])
        rr = HaltonN(Nm_slice, iprime = seed)
        
        current1 = 0
        Ne_slice = current1*lambda0/g_c/g_qe
        
        theta1 = theta1 = (2*rr-1)*np.pi
        
        Nm_slices += Nm_slice
        
        #print('# of slices: ', islice, ', current = ', current1, ' A')
        
        xx = temp0[:,0]*ccc0
        yy = temp0[:,2]*ccc0
        xp = temp0[:,1]*ccc0
        yp = temp0[:,3]*ccc0
        #theta1 = temp[:,4]-x1
        gamma = temp0[:,-1]
        
        # bf2 = np.abs(np.mean(np.exp(1j*theta1)))**2
        # stat.append([islice, x1, current1, bf2, bf2*(current1)**2, Ne_slice])
        
        bf = np.mean(np.exp(1j*theta1))
        bf2 = np.abs(bf)**2
        phase = np.angle(bf)
        stat.append([islice, x1, current1, bf2, bf2*(current1)**2, Ne_slice, 
                     phase, np.real(bf), np.imag(bf)])
        
        slicename = "/slice%06d/" % (islice+1)
        f.create_dataset(slicename+"current", data = [0])
        f.create_dataset(slicename+"gamma", data = gamma)
        f.create_dataset(slicename+"theta", data = theta1)
        f.create_dataset(slicename+"x", data = xx)
        f.create_dataset(slicename+"y", data = yy)
        f.create_dataset(slicename+"px", data = xp)
        f.create_dataset(slicename+"py", data = yp)
        
        
    stat = np.array(stat); curpeak = np.max(stat[:,2])
    header = str.format('%12s%14s%14s%14s' % ('islice', 'z', 'I', 'bf2'))
    with open('bunching_%d@' % Nm_slice+outputName, 'wb') as f_handle:
        np.savetxt(f_handle, stat, fmt = '%14.6E', 
                   header = header)
    
    f.create_dataset('/slicespacing', data = [lambda0])
    f.create_dataset('/slicelength', data = [lambda0])
    f.create_dataset('/refposition', data = [0])
    f.create_dataset('/slicecount', data = [Ns])
    nbins = 16
    
    f.create_dataset('/beamletsize', data = [nbins])
    if version == 4:
        f.create_dataset('/one4one', data = [0])
    
    f.close()
    
    ilast += 1
    fig, ax = plt.subplots(nrows = 4, figsize = (5, 9), sharex = True)
    ax1, ax2, ax3, ax4 = ax.flatten()
    ax1.plot(stat[:ilast,0], stat[:ilast,3], '-')
    ax1.plot(stat[:ilast,0], 1/stat[:ilast,5], 'k--')
    ax1.set_yscale('log')
    ax1.grid()
    #ax1.set_xlabel(r'Number of slices')
    ax1.set_ylabel(r'$|b_f|^2$')
    ax1.set_ylim(1e-10, 1)
    
    # ax2.plot(stat[:ilast,0], stat[:ilast,-3], '-'); 
    # print('Sum of bunching: ', np.sum(stat[:ilast,-3]))
    # ax2.set_yscale('log')
    # ax2.set_ylabel(r'$N^2|b_f|^2$')
    
    #ax2.plot(stat[:ilast,0], stat[:ilast,8]/stat[:ilast,7], '-')
    ax2.plot(stat[:ilast,0], stat[:ilast,6], '-')
    ax2.set_ylabel('Phase (rad)')
    
    #ax3.plot(stat[:ilast,0], stat[:ilast,7], '-')
    #ax3.plot(stat[:ilast,0], stat[:ilast,8], '-')
    
    ax3.plot(stat[:ilast,0], stat[:ilast,3]*stat[:ilast,2], '-')
    ax3.set_ylabel(r'$I|b_f|^2$')
    ax3.set_yscale('log')
    
    #ax4.plot(stat[:ilast,1]*1e3, stat[:ilast,2], '-')
    ax4.plot(stat[:ilast,0], stat[:ilast,2], '-')
    ax4.grid()
    ax4.set_xlabel(r'$z$ (mm)')
    ax4.set_ylabel(r'$I_{peak}$ (A)')
    ax4.set_xlabel(r'Number of slices')
    
    fig.tight_layout()
    fig.savefig('bunching@'+outputName+"_nperlambda-"+str(nperlambda)+'.png')
    
    return curpeak, curlen, Ns, nbins, outputName, beam_P0

def Astra2Genesis3Slices(inputName = None, inputDist = None, 
                         outputName = 'temp',
                         seed = 8, npart = 8192*4, nslice = None, 
                         lambda0 = 100e-6, degree = 4, nperlambda = 1,
                         Qscale = 1, zscale = 1, bunch = 0,
                         useHammersley = 1):
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
    bunch : float, optional
        Bunching factor introduced evenly to each slice. The default is zero.
    useHammersley: int, optional
        If 1, random samples will be generated using Hammersley sequence, otherwise 
        using the method proposed by C. Penman, 1992. The default is 1.
        
    Returns
    -------
    Ntail : int
        Same as `Ntail` in `Genesis1.3`.
    curpeak : float
        Peak current, A.
    curlen : float
        RMS length of the current profile, positive for Gaussian and negative for flattop,
        meter.
    Ns: int
        Same as `nslice` in Genesis1.3 V2
    Ntail: int
        Same as `ntail` in Genesis1.3 V2
    outputName: string
        Name of the file where all slices are saved, to be used as `partfile` in Genesis1.3 V2

    '''
    curpeak, curlen, Ns, nbins, outputName =  Astra2GenesisSlices(
        inputName = inputName, inputDist = inputDist, 
        outputName = outputName,
        seed = seed, npart = npart, nslice = nslice, 
        lambda0 = lambda0,
        Qscale = Qscale, zscale = zscale, bunch = bunch,
        useHammersley = useHammersley,
        version = 3)
    return curpeak, curlen, Ns, nbins, outputName

Astra2Genesis4Slices = Astra2GenesisSlices

#%% Test script of tapering undulator
def TaperedUndulator0(lam_u = 3e-2, Nu = 113, dK_dz = 0, d2K_dz2 = 0):
    
    Ku = 3.49 #
    
    Nu0 = 0
    
    ss = '!LOOP=2'
    print(ss)
    
    ss = str.format('%-6s%12.4E%12.4E%6d' % ('AW', 0, lam_u*3.5, 1))
    print(ss)
    
    Dz = (Nu-Nu0)*lam_u
    DKu = dK_dz*Dz+d2K_dz2*Dz**2/2
    #Ku0 = Ku-DKu
    Ku0 = Ku
    
    KK = []
    for i in np.arange(Nu-Nu0):
        zi = i*lam_u
        Ki = (Ku0+dK_dz*zi+d2K_dz2*zi**2/2)/np.sqrt(2)
        
        ss = str.format('%-6s%12.4E%12.4E%6d' % ('AW', Ki, lam_u, 1))
        print(ss)
    
        KK.append([zi, Ki])
    
    for i in np.arange(Nu0):
        zi = lam_u*(Nu-Nu0+i)
        Ki = Ku/np.sqrt(2)
        
        ss = str.format('%-6s%12.4E%12.4E%6d' % ('AW', Ki, lam_u, 1))
        print(ss)
        
        KK.append([zi, Ki])
    
    KK = np.array(KK)
    fig, ax = plt.subplots()
    ax.plot(KK[:,0], KK[:,1]/np.max(KK[:,1]))
    
    ss = str.format('%-6s%12.4E%12.4E%6d' % ('AW', 0, lam_u*3.5, 1))
    print(ss)
    
    ss = '!ENDLOOP'
    print(ss)

def TaperedUndulator(lam_u = 3e-2, Nu = 113, dK_dz = 0, d2K_dz2 = 0):
    
    Ku = 3.49 #
    
    Nu0 = 72-3
    Nu1 = 90-3
    
    ss = '!LOOP=2'
    print(ss)
    
    ss = str.format('%-6s%12.4E%12.4E%6d' % ('AW', 0, lam_u*3.5, 1))
    print(ss)
    
    #Dz = (Nu1-Nu0+1)*lam_u
    #DKu = dK_dz*Dz+d2K_dz2*Dz**2/2
    #Ku0 = Ku-DKu
    Ku0 = Ku
    
    aa = []
    
    for i in np.arange(Nu0):
        zi = i*lam_u
        ai = Ku0/np.sqrt(2)
        
        ss = str.format('%-6s%12.4E%12.4E%6d' % ('AW', ai, lam_u, 1))
        print(ss)
    
        aa.append([zi, ai])
        
    for i in np.arange(Nu0, Nu1):
        zi = (i-Nu0)*lam_u
        ai = (Ku0+dK_dz*zi+d2K_dz2*zi**2/2)/np.sqrt(2)
        
        ss = str.format('%-6s%12.4E%12.4E%6d' % ('AW', ai, lam_u, 1))
        print(ss)
    
        aa.append([zi+Nu0*lam_u, ai])
    
    au1 = ai; 
    for i in np.arange(Nu1, Nu):
        zi = lam_u*i
        Ki = au1
        
        ss = str.format('%-6s%12.4E%12.4E%6d' % ('AW', Ki, lam_u, 1))
        print(ss)
        
        aa.append([zi, Ki])
    
    aa = np.array(aa)
    fig, ax = plt.subplots()
    ax.plot(aa[:,0], aa[:,1]/np.max(aa[:,1]))
    
    ss = str.format('%-6s%12.4E%12.4E%6d' % ('AW', 0, lam_u*3.5, 1))
    print(ss)
    
    ss = '!ENDLOOP'
    print(ss)
#TaperedUndulator(dK_dz = -0.1)
        
#%% Main function
if __name__ == "__main__":
    import os
    # workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim3\2023\THzApple2\P0_15MeV_c_lam_115um_betax_0.3m'
    # workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim3\2023\SASEScan2'
    # workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim3\2023\THzDLW'
    # workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim3\2023\THzComm2'
    # workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim3\2024\WaveGuide'
    #workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim3\2024\S2E-Georgi'
    
    #workdir = r'\\afs\ifh.de\group\pitz\data\lixiangk\sim3\2024\DLW-XYZ'
    workdir = r'/mnt/d/NextCloud/subjects/2025/202501-2nC/debug_gen4_bugs/bug_version'   
    os.chdir(workdir)
    
    #fname = 'beam.ini'
    #fname = '368A.2809.002.1000'
    #fname = 'ast.2745.019_14.7_0.3'
    #fname = 'ast_1.0nC.2809.1500' 
    #fname = '.'+os.sep+'ast.2745.013'
    #fname = '22MeV_c.input.1.0'
    #fname = 'ast.2809.002.matched2'
    #fname = 'ast.2809.002.matched2.compressed'
    #fname = 'ast.2745.019'
    #fname = 'ast.2809.002'
    #fname = 'ast.2809.002.matched2_.compressed'
    #fname = 'g_100A_2.0nC.ini'
    #fname = 'run1.2809.001'
    #fname = 'DLW-Chicane-Und.250K.2809.001'
    #fname = 'DLW-Und.250K.2809.001'
    
    #fname = 'oce_2.5nC.2809.001'
    
    #fname = 'beam_112A_2.0nC.ini'
    fname = 'fort.1032.001'
    
    #Ipeak, Qtot = 200, 1
    # curpeak, curlen, Nslice, Ntail, outputName = Astra2Genesis3Slices(fname, 
    #                                               outputName = 'test', 
    #                                               seed = 8, 
    #                                               npart = 4096*2,
    #                                               Qscale = 1, 
    #                                               zscale = 1,
    #                                               bunch = 0)
    
    curpeak, curlen, Nslice, Nbins, outputName, beam_P0 = Astra2GenesisSlices(fname, 
                                                  outputName = 'scan', 
                                                  seed = 5, 
                                                  npart = 4096,
                                                  Qscale = 1, #Qtot/2, 
                                                  zscale = 1, #1.115, #(Qtot/Ipeak)/(2/100),
                                                  lambda0 = 100e-6,
                                                  degree = 1,
                                                  nperlambda = 3,
                                                  bunch = 0.0, 
                                                  useHammersley = 0)
    # outputName = Astra2WarpHaltonZ(fname, 
    #                                outputName = 'temp', 
    #                                seed = 5, 
    #                                Qscale = 1, #Qtot/2, 
    #                                zscale = 1, #(Qtot/Ipeak)/(2/100),
    #                                lambda0 = 100e-6,
    #                                bunch = 0.0, 
    #                                useHammersley = 1)
    
    ###

'''
from universal import easy_smooth
            def empirical_cdf_inv(samples):
                sorted_samples = np.sort(samples)
                n = len(samples)
                cdf_values = np.arange(1, n + 1) / n
                
                sorted_samples = easy_smooth(sorted_samples, 17)
                f_cdf_inv = sp.interpolate.interp1d(cdf_values, sorted_samples,
                                                    bounds_error = False,
                                                    fill_value = 'extrapolate')
                return f_cdf_inv
            
            
            def empirical_cdf_inv_(samples):
                sorted_samples = np.sort(samples)
                n = len(samples)
                cdf_values = np.arange(1, n + 1) / n
                
                # f_cdf_inv = sp.interpolate.interp1d(cdf_values, sorted_samples,
                #                                     bounds_error = False,
                #                                     fill_value = 'extrapolate')
                
                bins = 64
                weights, edges = np.histogram(samples, bins = bins)
                ### Add zeros to avoid abrupt change of density near tails
                #weights = np.concatenate(([0]*3, weights, [0]*3))
                dz = edges[1]-edges[0]
                #edges = np.concatenate((np.array([-3, -2, -1])*dz+edges[0], 
                #                        edges, np.array([1, 2, 3])*dz+edges[-1]))
                
                weights = easy_smooth(weights, 7)
                
                edges, weights = edges[1:], weights
                
                weights[0] = 0
                for j in np.arange(1, len(weights)):
                    weights[j] += weights[j-1]
                #weights = np.concatenate(([0], weights))
                weights = weights/weights[-1]
                ###

    
    
                Fx_inv = interp1d(weights, edges, kind = 'cubic')
                
                return Fx_inv
'''
