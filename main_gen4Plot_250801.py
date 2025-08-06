import h5py
import xkinterface.startup as stp
import os
import matplotlib.pyplot as plt
import numpy as np

stp.setplot()

try:
    hid.close()
except:
    pass

path = stp.winpath(r"F:\simu_2025\initialSamplingMethod\plannerUndu\04_match2twiss_timeON_SmoothGaussProfile\test_timeOFF")
path="/mnt/f/simu_2025/gen4_examples/Example1-SteadyState"
path="/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/04_match2twiss_timeON_SmoothGaussProfile/test_timeOFF"

path="/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/04_match2twiss_timeON_SmoothGaussProfile/test_timeOFF_taperUndu"

path="/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/04_match2twiss_timeON_SmoothGaussProfile/test_timeOFF_FODO"

path='/mnt/c/Users/biaob/Desktop/02_addFODO/test_timeOFF_FODO'

path='/mnt/c/Users/biaob/Desktop/02_addFODO/test_timeOFF_xyFocus'

path=stp.winpath(r"F:\simu_2025\initialSamplingMethod\plannerUndu\04_match2twiss_timeON_SmoothGaussProfile")
path=stp.winpath(r"F:\simu_2025\initialSamplingMethod\plannerUndu\04_xyFocus_smoothGaussProfile\test_timeON_xyFocus")

# path=stp.winpath(r"F:\simu_2025\initialSamplingMethod\plannerUndu\04_match2twiss_timeON_SmoothGaussProfile\SC_ON")

# path="/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/04_match2twiss_timeON_SmoothGaussProfile/test_timeOFF_2"

path="/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/06_one4one_smoothFlattopProfile/01_gen4"
path="/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/06_one4one_smoothFlattopProfile/01_gen4_gainCurveMeasure"

path=stp.winpath(r"F:\simu_2025\initialSamplingMethod\plannerUndu\04_match2twiss_timeON_shotON")
# path=stp.winpath(r"F:\simu_2025\initialSamplingMethod\plannerUndu\04_match2twiss_timeON_shotON\test_radiationProfileEvo")

# path="/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/06_one4one_smoothFlattopProfile/01_gen4_shotNoiseON"
os.chdir(path)

# fname = 'Example1.out.h5'
fname = "gen4.out.h5"
hid = h5py.File(fname,'r')

#%% beam info, energy spread and detune
from xkinterface.interface import postgen4

# postgen4.plot_BeamEnergyDetune(hid)
# stp.savefig()

# postgen4.plot_BeamEnergySpread(hid)
# stp.savefig()


# z = hid['Lattice']['z'][::3]
bE0 = hid["Beam"]["Global"]["energy"][:] *0.511 #MeV
sigE = hid["Beam"]["Global"]["energyspread"][:] *0.511  #MeV

tmpN = len(bE0)
sss = np.linspace(0, 3.6, tmpN)


E0_i = bE0[0]
eta_off = (bE0-E0_i)/E0_i

plt.figure()
# plt.plot(sss,eta_off*100,'-')
plt.plot(sss, bE0,'-')

plt.xlabel("s (m)")
plt.ylabel("beam energy (MeV)")
# plt.ylabel("beam energy detune (%)")
plt.grid()

stp.savefig("beam_detune")

#%%% energy change of given slice
slil = np.arange(0,100,5)
zz = hid["Lattice"]["zplot"][:]

plt.figure()
for sli in slil:
    
    sli_E0 = hid["Beam"]["energy"][:,sli] *0.511 #MeV
    
    plt.plot(zz, sli_E0, '-',label=f"slice {sli}")

plt.legend()    
plt.grid()
stp.savefig("slice_energy")
    
#%% beam current
from xkinterface.interface import postgen4

current = postgen4.plot_BeamCurrentProfile(hid)
stp.savefig()


plt.figure()
plt.plot(current,'-')
plt.xlabel("slice #")
plt.ylabel("current (A)")
plt.grid()

stp.savefig("slice_current")

#%% local bunching
jj=0
spos = 0.015*3*jj


bf_local = hid["Beam"]["bunching"][jj,:] 

plt.figure()

plt.title(f"s={spos:.1f} m")
s = hid['Global']['s'][()]/2.998e8*1e12  #/5.5 -2.97 #-4.9 #ps

plt.plot(s, bf_local,'--')
# plt.xlim([0,30])
# plt.ylim([1e-5,1])

plt.semilogy()
plt.grid()

#%%% local bunching evolution
jjl = np.arange(0,75,10)

plt.figure()
for jj in jjl:
    spos = 0.015*jj*3

    bf_local = hid["Beam"]["bunching"][jj,:] 
    s = hid['Global']['s'][()]/2.998e8*1e12  #/5.5 -2.97 #-4.9 #ps
    
    plt.plot(s, bf_local,'-',label=f"s={spos:.1f} m")
    # plt.xlim([-3,3])
    # plt.ylim([1e-5,1])

    plt.semilogy()
    plt.grid()
plt.legend()
plt.xlabel("z (ps)")
plt.ylabel("local bunching")

stp.savefig()

#%% gain curve
from xkinterface.interface import postgen4

postgen4.plot_gainCurvePeakPower(hid,semilog=True)
stp.savefig(filename="gainCurve_peakPower")

postgen4.plot_gainCurvePower(hid,semilog=True)
stp.savefig(filename="gainCurve_rmsSize")


#%% rms size
from xkinterface.interface import postgen4

postgen4.plot_rmsSizeEvo(hid, field=False)
stp.savefig(filename="rmsSize")

postgen4.plot_rmsSizeEvo(hid, field=True)
stp.savefig(filename="rmsSize_field")

#%% power profile evo
from xkinterface.interface import postgen4

postgen4.plot_powerEvolution(hid, norm=False)
stp.savefig()

postgen4.plot_powerEvolution(hid, norm=True)
stp.savefig("powerProfile")

#%% beam slice energy evo
from xkinterface.interface import postgen4

postgen4.plot_beamEnergyEvolution(hid, sliceRange=[0,90])
stp.savefig("sliceEnergy")

#%% beam slice bunching evo
from xkinterface.interface import postgen4

postgen4.plot_beamBunchingEvolution(hid, norm=True, sliceRange=[5,90])
stp.savefig("sliceBunchingNorm")

postgen4.plot_beamBunchingEvolution(hid, norm=False, sliceRange=[5,90])
stp.savefig("sliceBunching")
#%% spectrum
from xkinterface.interface import postgen4

# postgen4.plot_spectrum(hid, step=70)
# stp.savefig()

postgen4.plot_spectrumEvolution(hid, norm=False, freqRange=[2.5,3.5])
stp.savefig("spectrumEvo")

postgen4.plot_spectrumEvolution(hid, norm=True, freqRange=[2.5,3.5])
stp.savefig("spectrumEvo_norm")

#%%% spectrum heat plot



