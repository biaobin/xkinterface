import h5py
import numpy as np
import matplotlib.pyplot as plt
from xkinterface.startup import *

setplot()

def getBunchingFactor(hid, plot=True, fig=True, figext=False) :

    z = hid['Lattice']['zplot'][()]
    current = hid["Beam"]["current"]
    a = hid["Beam"]["bunchingreal"]
    b = hid["Beam"]["bunchingimag"]
    
    steps = a.shape[0]
    
    bunchingFac = []
    for j in range(steps):
        tmp1 = a[j,:] * current[0,:]
        tmp2 = b[j,:] * current[0,:]
        abs_bf = np.sqrt(np.sum(tmp1)**2 + np.sum(tmp2)**2)
        tmp = abs_bf/np.sum(current[0,:])

        bunchingFac.append( tmp )
    bunchingFac = np.array(bunchingFac)

    if plot==True:
        if fig==True:
            fig, ax1 = plt.subplots(figsize=(8,6))
            
        ax1.plot(z, bunchingFac,'-r', label='bunching factor')
        plt.grid()
        plt.xlabel('s (m)')
        plt.ylabel('bunching factor')
    
        if figext==True:
            savefig(filename="bunchingFac")
    
    return z,bunchingFac


def plot_gainCurvePeakPower(hid, plot=True, fig=True, line='-', figext=False, semilog=True):
    
    z = hid["Lattice"]["zplot"][:]
    peak_power = np.max( hid["Field"]["power"], axis=1)
    
    z, bunchingFac = getBunchingFactor(hid, plot=False)

    if plot==True:
        if fig==True:
            fig,ax1 = plt.subplots()
        
        if semilog == True:
            ax1.semilogy(z,peak_power/1e6,line)
        else:
            ax1.plot(z,peak_power/1e6,line)
        ax1.set_xlabel(r'$s$ (m)')
        ax1.set_ylabel(r'peak power (MW)')
        plt.grid(True)
        
        ax2 = ax1.twinx()
        color = 'tab:blue'
        ax2.set_ylabel(r'bunching factor',color = color)
        ax2.tick_params(axis='y', labelcolor=color)
        if semilog == True:
            ax2.semilogy(z,bunchingFac,line,color = color, label='bunching factor')
        else:
            ax2.plot(z,bunchingFac,line,color = color, label='bunching factor')
        plt.show()
        
        if figext==True:
            savefig(filename="gaincurvePeakPower")        

    return z, peak_power
    

def plot_gainCurvePower(hid, bunching=True, semilog=True, figext=False):

    z, bunchingFac = getBunchingFactor(hid, plot=False)
    
    z = hid['Lattice']['zplot'][()]
    energy = hid['Field']['Global']['energy'][()]*1e6  #uJ
    
    #nearfield = hid['Field']['Global']['intensity-nearfield'][()]*1e-21  #inital unit [w/rad]
    #farfield  = hid['Field']['Global']['intensity-farfield'][()]*1e-21  #inital unit [w/rad]  => GW/urad^2
    
    fig,ax1 = plt.subplots()
    color = 'tab:red'
    ax1.set_xlabel(r'$z$ (m)')
    ax1.set_ylabel(r'$E$ ($\mu$J)',color = color)
    ax1.tick_params(axis='y', labelcolor=color)
    if semilog == True:
        ax1.semilogy(z,energy,'-',color = color,label='energy')
    else:
        ax1.plot(z,energy,'-',color = color,label='energy')
    # ax1.set_ylim([1e-3,1e3])
    # plt.legend()
    plt.grid(True)
    
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel(r'bunching factor',color = color)
    ax2.tick_params(axis='y', labelcolor=color)
    if semilog == True:
        ax2.semilogy(z,bunchingFac,'-',color = color, label='bunching factor')
    else:
        ax2.plot(z,bunchingFac,'-',color = color, label='bunching factor')
    # ax2.semilogy(z,farfield,color ='g', label='farfield')
    # ax2.set_ylim([1e-2,1e5])
    # plt.legend()
    plt.show()
    
    if figext==True:
        savefig(filename="gaincurve")

def plot_powerEvolution(hid, norm=True, figext=False):
    
    z = hid['Lattice']['zplot'][:]
    s = hid['Global']['s'][()]/2.998e8*1e12  #ps
    power =  hid['Field']['power'][()]
    
    fig,ax1 = plt.subplots()
    if norm == True:
        plt.title(r"$power_{step_i}/mean(power_{step_i})$")
    
        pmean= np.mean(power, axis = 1)
        for i in range(len(pmean)):
            if pmean[i] == 0:
                pmean[i]=1.
            power[i,:]*=1./pmean[i]
        plt.imshow(np.flipud(power), aspect='auto', interpolation='none', extent=(np.min(s),np.max(s),np.min(z),np.max(z)))
        plt.xlabel(r'z (ps)')
        plt.ylabel(r'$s$ (m)')
        # plt.colorbar()
        plt.show()
    else:
        plt.title("power evolution (W)")
        # plt.imshow(np.flipud(np.log10(abs(power))), aspect='auto', interpolation='none',extent=(np.min( s ),np.max(s),np.min( z ),np.max(z)))
        plt.imshow(np.flipud(power), aspect='auto', interpolation='none',extent=(np.min( s ),np.max(s),np.min( z ),np.max(z)))
        plt.xlabel(r'z (ps)')
        plt.ylabel(r'$s$ (m)')
        plt.colorbar()
        plt.show()
    
    if figext==True:
        savefig(filename="powerevo")
        
        
def plot_beamEnergyEvolution(hid, figext=False, sliceRange=[0,90]):
    
    z = hid['Lattice']['zplot'][:]
    
    #beam energy heat plot
    s = hid['Global']['s'][sliceRange[0]:sliceRange[1]]/2.998e8*1e12  #ps
    power = hid["Beam"]["energy"][:,sliceRange[0]:sliceRange[1]]*0.511 #MeV

    fig,ax1 = plt.subplots()
    plt.title("beam slice energy evolution (MeV)")
    # plt.imshow(np.flipud(np.log10(abs(power))), aspect='auto', interpolation='none',extent=(np.min( s ),np.max(s),np.min( z ),np.max(z)))
    plt.imshow(np.flipud(power), aspect='auto', interpolation='none',extent=(np.min( s ),np.max(s),np.min( z ),np.max(z)))
    plt.xlabel(r'z (ps)')
    plt.ylabel(r'$s$ (m)')
    plt.colorbar()
    plt.show()
    
    if figext==True:
        savefig(filename="beamSliceEnergyEvo")
        

def plot_beamBunchingEvolution(hid, figext=False, norm=False, sliceRange=[0,90]):
    
    z = hid['Lattice']['zplot'][:]
    
    #beam energy heat plot
    s = hid['Global']['s'][sliceRange[0]:sliceRange[1]]/2.998e8*1e12  #ps
    power = hid["Beam"]["bunching"][:,sliceRange[0]:sliceRange[1]]

    if norm == True:
        pmean= np.mean(power, axis = 1)
        for i in range(len(pmean)):
            if pmean[i] == 0:
                pmean[i]=1.
            power[i,:]*=1./pmean[i]

    fig,ax1 = plt.subplots()
    if norm == True:
        plt.title("beam slice norm. bunching evolution")
    else:
        plt.title("beam slice bunching evolution")
    # plt.imshow(np.flipud(np.log10(abs(power))), aspect='auto', interpolation='none',extent=(np.min( s ),np.max(s),np.min( z ),np.max(z)))
    plt.imshow(np.flipud(power), aspect='auto', interpolation='none',extent=(np.min( s ),np.max(s),np.min( z ),np.max(z)))
    plt.xlabel(r'z (ps)')
    plt.ylabel(r'$s$ (m)')
    plt.colorbar()
    plt.show()
    
    if figext==True:
        savefig(filename="beamSliceEnergyEvo")
        
def plot_spectrumEvolution(hid, norm=False, figext=False, freqRange=None):
    energy = hid['Field']['Global']['energy'][()]*1e6
    steps = len(energy)
        
    specll = []
    for step in range(steps):
        
        # print(step)
        energy0=energy[step]
        
        sig = hid['Field']['intensity-farfield'][()][step,:]
        phi = hid['Field']['phase-farfield'][()][step,:] 
        freq = hid['Global']['frequency'][()]
        
        signal = np.sqrt(sig)*np.exp(1j*phi)
        spec = np.abs(np.fft.fftshift(np.fft.fft(signal)))**2
        norm_v = energy0/np.sum(spec)/(freq[1]-freq[0])
        spec =norm_v*spec
        
        specll.append(spec)
        # plt.plot(f_THz, spec,'-')
    
    specll = np.array(specll)
    
    if norm == True:
        pmean= np.mean(specll, axis = 1)
        for i in range(len(pmean)):
            if pmean[i] == 0:
                pmean[i]=1.
            specll[i,:]*=1./pmean[i]
    
    z = hid['Lattice']['zplot'][:]
    # frequency, eV => Hz
    # 1eV = 1.6e-19 J, h=6.626e-34
    f_THz = freq*1.602e-19/6.626e-34/1e12        
    
    if freqRange == None:
        freqRange = [np.min(f_THz), np.max(f_THz)]
    else:
        
        idd = (f_THz > freqRange[0]) & (f_THz < freqRange[1])
        
        f_THz = f_THz[idd]
        specll = specll[:,idd]
    
    fig,ax1 = plt.subplots()
    
    # plt.imshow(np.flipud(np.log10(abs(specll))), aspect='auto', interpolation='none',extent=(np.min( s ),np.max(s),np.min( z ),np.max(z)))
    plt.imshow(np.flipud(specll), aspect='auto', interpolation='none',extent=(np.min(f_THz),np.max(f_THz),np.min(z),np.max(z)))
    plt.xlabel(r'frequency (THz)')
    plt.ylabel(r'$s$ (m)')
    
    if norm == False:
        plt.colorbar()
        plt.title("beam spectrum evolution")
    else:
        plt.title("norm. beam spectrum evolution")
    plt.show()
    
    if figext==True:
        savefig(filename="spectrumEvo")

def plot_spectrum(hid, step=-1, fig=True):
    
    energy = hid['Field']['Global']['energy'][()]*1e6
    energy0=energy[step]
    
    sig = hid['Field']['intensity-farfield'][()][step,:]
    phi = hid['Field']['phase-farfield'][()][step,:] 
    freq = hid['Global']['frequency'][()]

    signal = np.sqrt(sig)*np.exp(1j*phi)
    spec = np.abs(np.fft.fftshift(np.fft.fft(signal)))**2
    norm = energy0/np.sum(spec)/(freq[1]-freq[0])
    spec =norm*spec

    # frequency, eV => Hz
    # 1eV = 1.6e-19 J, h=6.626e-34
    f_THz = freq*1.602e-19/6.626e-34/1e12

    # plt.figure()

    # plt.title(f"s={spos:.2f} m")
    # plt.plot(freq*1e-3,spec)
    # # plt.xlim([12.400,12.800])
    # plt.xlabel(r'$E_{ph}$ (keV)')
    # plt.ylabel(r'$P(E_{ph})$ ($\mu$J/eV)')
    # plt.show()
    
    if fig==True:
        plt.figure()
    # plt.title(f"s={spos:.2f} m")
    
    plt.plot(f_THz,spec,'-')
    # plt.xlim([12.400,12.800])
    plt.xlabel(r'$f_{ph}$ (THz)')
    plt.ylabel(r'$P(E_{ph})$ ($\mu$J/eV)')
    plt.show()

        
def plot_rmsSizeEvo(hid, output_step=1, fig=True, field=True):

    z = hid['Lattice']['zplot'][::output_step]
    bx = hid['Beam']["Global"]['xsize'][:]
    by = hid['Beam']["Global"]['ysize'][:]
    fx = hid['Field']["Global"]['xsize'][:]
    fy = hid['Field']["Global"]['ysize'][:]

    if fig == True:
        plt.figure()
    plt.plot(z,bx*1e3,'-',label=r'Beam: $\sigma_x$')
    plt.plot(z,by*1e3,'-',label=r'Beam: $\sigma_y$')
    if field==True:
        plt.plot(z,fx*1e3,'-',label=r'Field: $\sigma_x$')
        plt.plot(z,fy*1e3,'-',label=r'Field: $\sigma_y$')
    plt.legend()
    plt.xlabel(r'$s$ (m)')
    plt.ylabel(r'$\sigma_{x,y}$ (mm)')
    # plt.ylim([0,60])

    plt.grid()
    
def plot_rmsSizeEvo_sliceJ(hid, slicej=10 ,output_step=1, fig=True, field=True):
    
    if fig==True:
        plt.figure()

    z = hid['Lattice']['zplot'][()]
    bx = hid['Beam']['xsize'][:,slicej]
    by = hid['Beam']['ysize'][:,slicej]
    fx = hid['Field']['xsize'][:,slicej]
    fy = hid['Field']['ysize'][:,slicej]

    plt.plot(z,bx*1e3,'-',label=r'Beam: $\sigma_x$')
    plt.plot(z,by*1e3,'-',label=r'Beam: $\sigma_y$')
    
    if field==True:
        plt.plot(z,fx*1e3,'-',label=r'Field: $\sigma_x$')
        plt.plot(z,fy*1e3,'-',label=r'Field: $\sigma_y$')
    plt.legend()
    plt.xlabel(r'$s$ (m)')
    plt.ylabel(r'$\sigma_{x,y}$ (mm)')
    # plt.ylim([0,60])

    plt.grid()

    plt.show()
    

def plot_BeamEnergyDetune(hid, sampleFreq=3, fig=True):
    # z = hid['Lattice']['z'][::sampleFreq]
    # bE0 = hid["Beam"]["Global"]["energy"][0:-1] *0.511 #MeV
    # sigE = hid["Beam"]["Global"]["energyspread"][0:-1] *0.511  #MeV
    
    z = hid['Lattice']['z'][:]
    bE0 = hid["Beam"]["Global"]["energy"][:] *0.511 #MeV
    sigE = hid["Beam"]["Global"]["energyspread"][:] *0.511  #MeV

    E0_i = bE0[0]
    eta_off = (bE0-E0_i)/E0_i

    if fig==True:
        plt.figure()

    plt.plot(z,eta_off*100,'-')

    plt.xlabel("s (m)")
    plt.ylabel("beam energy detune (%)")
    plt.grid()


def plot_BeamEnergySpread(hid, sampleFreq=3, fig=True):
    # z = hid['Lattice']['z'][::sampleFreq]
    # bE0 = hid["Beam"]["Global"]["energy"][0:-1] *0.511 #MeV
    # sigE = hid["Beam"]["Global"]["energyspread"][0:-1] *0.511  #MeV
    
    z = hid['Lattice']['z'][:]
    bE0 = hid["Beam"]["Global"]["energy"][:] *0.511 #MeV
    sigE = hid["Beam"]["Global"]["energyspread"][:] *0.511  #MeV

    E0_i = bE0[0]
    eta_off = (bE0-E0_i)/E0_i

    if fig==True:
        plt.figure()
        
    plt.plot(z,sigE/bE0*100,'-')
    plt.xlabel("s (m)")
    plt.ylabel(r"beam energy spread $\Delta E/E$ (%)")
    plt.grid()


def plot_BeamCurrentProfile(hid):
    s = hid['Global']['s'][()]/2.998e8*1e12 #ps
    current = hid['Beam']['current'][()][0,:]
    energy = hid['Beam']['energy'][()][0,:]*0.511


    fig,ax1 = plt.subplots()
    color = 'tab:red'
    ax1.set_xlabel(r'$s$ (ps)')
    ax1.set_ylabel(r'$I$ (A)',color = color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.plot(s,current,'-',color = color)
    plt.grid()

    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel(r'$E$ (MeV)',color = color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.plot(s,energy,'-',color = color)
    plt.show()
    
    return current 
    
    
    
    
