import h5py
import xkinterface.bbl.startup as stp
import os
import matplotlib.pyplot as plt
import numpy as np

stp.setplot()

try:
    hid.close()
except:
    pass

#%%
path='/mnt/f/simu_2025/202509_sase1D_genesisSimu_THzBenchmark/02_smoothGaussProfile'
path='/mnt/f/simu_2025/202509_sase1D_genesisSimu_THzBenchmark/01_shotnoise_0.25mm_120A_gen4_seeding'

path='/mnt/f/simu_2025/202510_idealMachineTHzFEL/test_1kA/Ipeak120A'
# path='/mnt/f/simu_2025/202510_idealMachineTHzFEL/Ipeak120_static'

path='/mnt/f/simu_2025/202510_idealMachineTHzFEL/05_bunchCompressorDesign/02_ocelot/02_gen4'
path='/mnt/f/simu_2025/202510_idealMachineTHzFEL/202511_scan_currentProfile_parabolic_flattop/00_ScanPythonScripts/Planar_flattop_slscON_llscON/Ipeak350A_Lbuncht6.40ps_Q2.02nC/02_gen4'
path='/mnt/f/simu_2025/202510_idealMachineTHzFEL/05_bunchCompressorDesign/202511_gaussComp/02_L2_dechirpScan/02_dechirped'
path='/mnt/f/simu_2025/202510_idealMachineTHzFEL/05_bunchCompressorDesign/202511_gaussComp/theta_11.0deg_300A'

path='/mnt/f/simu_2025/202510_idealMachineTHzFEL/05_bunchCompressorDesign/202511_gaussComp/02_L2_dechirpScan/03_dechirped_525A_colimYYP'

path='/mnt/f/simu_2025/202510_idealMachineTHzFEL/05_bunchCompressorDesign/02_ocelot/02_gen4/ave_y_yp_eq0/..'

path='/mnt/f/simu_2025/202510_idealMachineTHzFEL/0000_syncDESY/gauss_Q2nC_phiboo_-32deg/03_gen4_enx6_eny38_sigE60_Lbunch5.5ps_2nC/Planar_flattop_slscON_llscON_one4oneFalse_fixCharge/Q2.0nC_Lbuncht5.50ps/02_gen4'

path=stp.winpath(r'F:\simu_2025\202510_idealMachineTHzFEL\05_bunchCompressorDesign\202511_gaussComp\02_L2_dechirpScan\03_dechirped_525A_colimYYP')


path='/mnt/f/simu_2025/202510_idealMachineTHzFEL/0000_syncDESY/beamEnergy_17MeV/gauss_Q2nC_phiboo_-40deg/02_gen4/nocollim'
os.chdir(path)

#%%
# fname = 'Example4_b.out.h5'
fname = "gen4.out.h5"
# fname = 'pithz.5.out.h5'
# fname='g4.100.out.h5'

hid = h5py.File(fname,'r')

#%% space charge

tmp = hid["Beam"]["SSCfield"][10,:]


plt.figure()
plt.plot(tmp,'-')

#%% beam info, energy spread and detune
from xkinterface.interface import postgen4

postgen4.plot_BeamEnergyDetune(hid)
stp.savefig()

postgen4.plot_BeamEnergySpread(hid)
stp.savefig()


# # z = hid['Lattice']['z'][::3]
# bE0 = hid["Beam"]["Global"]["energy"][:] *0.511 #MeV
# sigE = hid["Beam"]["Global"]["energyspread"][:] *0.511  #MeV

# tmpN = len(bE0)
# zz = hid["Lattice"]["zplot"][:]
# sss = np.linspace(np.min(zz), np.max(zz), tmpN)


# E0_i = bE0[0]
# eta_off = (bE0-E0_i)/E0_i

# plt.figure()
# # plt.plot(sss,eta_off*100,'-')
# plt.plot(sss, bE0,'-')

# plt.xlabel("s (m)")
# plt.ylabel("beam energy (MeV)")
# # plt.ylabel("beam energy detune (%)")
# plt.grid()

# stp.savefig("beam_detune")

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
from xkinterface.bbl import postgen4

try:
    hid.close()
except:
    pass

fname = "gen4.out.h5"
hid = h5py.File(fname,'r')

current = postgen4.plot_BeamCurrentProfile(hid)
# stp.savefig()

current = postgen4.plot_BeamCurrentProfile(hid,slicex=True)
stp.savefig()

# plt.figure()
# plt.plot(current,'-')
# plt.xlabel("slice #")
# plt.ylabel("current (A)")
# plt.grid()

# stp.savefig("slice_current")

#%% local bunching
jj=0
spos = 0.1*jj


bf_local = hid["Beam"]["bunching"][jj,:] 

plt.figure()

plt.title(f"s={spos:.1f} m")
s = hid['Global']['s'][()]/2.998e8*1e12  #/5.5 -2.97 #-4.9 #ps

plt.plot(s, bf_local,'--')
# plt.xlim([0,30])
# plt.ylim([1e-5,1])

plt.semilogy()
plt.grid()


#%%local bunching evolution of given slice
steps_num = hid["Beam"]["bunching"].shape[0]
slice_num = hid["Beam"]["bunching"].shape[1]

slil = np.arange(0,slice_num,1)
zz = hid["Lattice"]["zplot"][:]

plt.figure()
for sli in slil:
    
    sli_E0 = hid["Beam"]["bunching"][:,sli]
    
    plt.plot(zz, sli_E0, '-',label=f"slice {sli}")

plt.legend()    
plt.grid()
plt.xlabel("s (m)")
plt.ylabel("slice bunching")

stp.savefig("slice_bunching")


#%%% local bunching evolution
steps_num = hid["Beam"]["bunching"].shape[0]
jjl = np.arange(0,steps_num,20)

plt.figure()
for jj in jjl:
    # spos = 0.015*jj*3
    spos = 0.015*jj*1

    bf_local = hid["Beam"]["bunching"][jj,:] 
    s = hid['Global']['s'][()]/2.998e8*1e12  #/5.5 -2.97 #-4.9 #ps
    
    # plt.plot(s[0:96], bf_local[0:96],'-',label=f"s={spos:.1f} m")
    plt.plot(s[0:200], bf_local[0:200],'-',label=f"s={spos:.1f} m")
    # plt.xlim([-3,3])
    # plt.ylim([1e-5,1])

    plt.semilogy()
    plt.grid()
plt.legend()
plt.xlabel("z (ps)")
plt.ylabel("local bunching")

stp.savefig()

#%% gain curve with fitting rho & Lg
# from xkinterface.interface import postgen4

# path='/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/04_match2twiss_timeON_shotON/gaussFromImpz/scOFF'

# path='/mnt/f/simu_2025/202509_XiangkunSimu/beam_112A_2.0nC_shotnoise_v4_quiet2'
# path='/mnt/f/simu_2025/202509_XiangkunSimu/beam_112A_2.0nC_shotnoise_v4_quiet2_2'

# path=stp.winpath(r'F:\simu_2025\initialSamplingMethod\plannerUndu\03_match2twiss_timeOFF_shotON')
# os.chdir(path)

try:
    hid.close()
except:
    pass

fname = "gen4.out.h5"
hid = h5py.File(fname,'r')

postgen4.plot_gainCurvePeakPower(hid,semilog=False)
stp.savefig(filename="gainCurve_peakPower")

postgen4.plot_gainCurvePower(hid,semilog=False)
stp.savefig(filename="gainCurve_energy")


# postgen4.plot_gainCurvePeakPower(hid,semilog=True)
# stp.savefig(filename="gainCurve_peakPower")

# postgen4.plot_gainCurvePower(hid,semilog=True)
# stp.savefig(filename="gainCurve_energy")

#%%% add fitting Lg and rho values
#-----------------------------
xx1 = hid["Lattice"]["zplot"][:]
# yy1 = np.max( hid["Field"]["power"], axis=1)
yy1 = hid['Field']['Global']['energy'][:]

x1,x2=0.9, 1.5
Lg, coeffs = postgen4.get_gainlen(xx1, yy1, zmin=x1, zmax=x2)

def get_rho(lamdau=30e-3, Lg=0.1):
    return lamdau/(4*np.pi*np.sqrt(3)*Lg)

rho = get_rho(Lg=Lg)

dx=0.0
xx = np.linspace(x1-dx, x2+dx, 50)
yy = np.exp(np.polyval(coeffs, xx))

plt.figure()
plt.plot(xx1/Lg, yy1 *1e6,'-')
plt.plot(xx/Lg,  yy*1e6,'--',label=fr'fitting, Lg={Lg:.4f} m, $\rho$={rho:.2e}')
plt.yscale('log')
plt.xlabel('s/Lg')
# plt.ylabel('peak power (MW)')
plt.ylabel('pulse energy (uJ)')

plt.legend()
plt.show()

stp.savefig(filename="gainCurve_peakPower_withLg")

#%%% bunching factor
z = hid['Lattice']['zplot'][()]
current = hid["Beam"]["current"]
a0 = hid["Beam"]["bunchingreal"]
b0 = hid["Beam"]["bunchingimag"]

A = np.array(hid["Beam"]["bunching"])
phi = np.array(hid["Beam"]["bunchingphase"])

a = np.sqrt(A**2/(1+np.tan(phi)**2))
a = np.sign(a0)*a
b = a*np.tan(phi)

steps = z.shape[0]

bunchingFac = []
for j in range(steps):
    tmp1 = a[j,:] * current[0,:]
    tmp2 = b[j,:] * current[0,:]
    abs_bf = np.sqrt(np.sum(tmp1)**2 + np.sum(tmp2)**2)
    tmp = abs_bf/np.sum(current[0,:])

    bunchingFac.append( tmp )
bunchingFac = np.array(bunchingFac)

# plt.figure()
plt.plot(z, bunchingFac,'--')

#%%% using XK's module
import xkinterface.interface.PostGenesis13 as PostGenesis
fname = "gen4.out.h5"

pg1 = PostGenesis(fname, debug=1, version=4, fig_ext = '.png')


#%%% compare
from xkinterface.interface import postgen4

# path=stp.winpath(r"F:\simu_2025\gen4_examples\Example4-HGHG\00_one4one_macroCharge1000")
# path = stp.winpath(r"F:\simu_2025\202508_HGHG_THz\pitzPara\02_HGHG_Np8192/scON")
# path='/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/04_match2twiss_timeON_shotON/sc_ON'

# os.chdir(path)

# fname = '/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/04_match2twiss_timeON_shotON/gaussFromImpz/scOFF/gen4.out.h5'
# fname=stp.winpath(r'F:\simu_2025\initialSamplingMethod\plannerUndu\04_match2twiss_timeON_shotON\scOFF')+'/gen4.out.h5'
# fname='/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/04_match2twiss_timeON_shotON/gaussFromImpz/scON_haltonImpz'+'/gen4.out.h5'

fname='/mnt/f/simu_2025/202509_XiangkunSimu/beam_112A_2.0nC_shotnoise_v4_quiet2'
fname += '/gen4.out.h5'

hid1 = h5py.File(fname,'r')

# fname = "/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/04_match2twiss_timeON_shotON/scOFF/gen4.out.h5"
# fname=stp.winpath(r'F:\simu_2025\initialSamplingMethod\plannerUndu\04_match2twiss_timeON_SmoothGaussProfile')+'/gen4.out.h5'
# fname='/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/04_match2twiss_timeON_shotON/scON_120A_150nwig'+'/gen4.out.h5'
fname='/mnt/f/simu_2025/202509_XiangkunSimu/beam_112A_2.0nC_shotnoise_v4_quiet2_'
fname += '/gen4.out.h5'

hid2 = h5py.File(fname,'r')

plt.figure()
postgen4.plot_gainCurvePeakPower2(hid1, semilog=True, figp=False, linestyle='-')
postgen4.plot_gainCurvePeakPower2(hid2, semilog=True, figp=False, linestyle='--')

plt.legend(['impz_gen','gen4_gen'])
plt.show()
stp.savefig(filename="compare_peakPower")

plt.figure()
postgen4.plot_gainCurvePower2(hid1, semilog=True, figp=False, linestyle='-')
postgen4.plot_gainCurvePower2(hid2, semilog=True, figp=False, linestyle='--')
# postgen4.plot_gainCurvePower2(hid3, semilog=True, figp=False, linestyle='--')

plt.legend(['impz_gen','gen4_gen'])
plt.show()
stp.savefig(filename="compare_energy")


# plt.figure()
# postgen4.plot_gainCurvePower2(hid, semilog=False, figp=False, linestyle='-')
# postgen4.plot_gainCurvePower2(hid2, semilog=False, figp=False, linestyle='--')
# stp.savefig(filename="compare_energy")


#%%% compare
from xkinterface.bbl import postgen4

# path = stp.winpath(r"F:\simu_2025\202508_HGHG_THz\pitzPara\02_HGHG_Np8192")
# path='/mnt/f/simu_2025/202510_idealMachineTHzFEL/test_1kA/Ipeak120A_scON'
# os.chdir(path)

# fname = '/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/04_match2twiss_timeON_shotON/scOFF/test/helical_beta_0.4_0.02/gen4.out.h5'
# hid = h5py.File(fname,'r')
# fname = "/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/04_match2twiss_timeON_shotON/scOFF/test/helical_beta_0.0922_0.0922193/gen4.out.h5"
# hid2 = h5py.File(fname,'r')

fname='/mnt/f/simu_2025/202510_idealMachineTHzFEL/0000_syncDESY/beamEnergy_17MeV/gauss_Q2nC_phiboo_-40deg/02_gen4/nocollim/gen4.out.h5'
hid = h5py.File(fname,'r')

fname='/mnt/f/simu_2025/202510_idealMachineTHzFEL/0000_syncDESY/beamEnergy_17MeV/gauss_Q2nC_phiboo_-40deg/02_gen4/collimY/gen4.out.h5'
hid2 = h5py.File(fname,'r')


fig,ax1 = plt.subplots()
postgen4.plot_gainCurvePeakPower2(hid, semilog=False, figp=False,linestyle='-',label='no alignment.')

postgen4.plot_gainCurvePeakPower2(hid2,semilog=False, figp=False,linestyle='--',label='alignment')
plt.legend()

stp.savefig('gain_compare')


plt.figure()
postgen4.plot_gainCurvePower2(hid, semilog=False, figp=False, linestyle='-',label='no alignment')
postgen4.plot_gainCurvePower2(hid2, semilog=False, figp=False, linestyle='--',label='alignment')
plt.legend()
plt.show()
stp.savefig(filename="compare_energy")


#%%% compare, chirp or non-dechirp
from xkinterface.bbl import postgen4

fname='/mnt/f/simu_2025/202510_idealMachineTHzFEL/0000_syncDESY/beamEnergy_40MeV/03_gen4/gen4.out.h5'
hid = h5py.File(fname,'r')

fname='/mnt/f/simu_2025/202510_idealMachineTHzFEL/0000_syncDESY/beamEnergy_40MeV/03_gen4_dechirp/gen4.out.h5'
hid2 = h5py.File(fname,'r')


fig,ax1 = plt.subplots()
postgen4.plot_gainCurvePeakPower2(hid, semilog=False, figp=False,linestyle='-',label='keep chirp')

postgen4.plot_gainCurvePeakPower2(hid2,semilog=False, figp=False,linestyle='--',label='rm chirp')
plt.legend()

stp.savefig('gain_compare')


plt.figure()
postgen4.plot_gainCurvePower2(hid, semilog=False, figp=False, linestyle='-',label='keep chirp')
postgen4.plot_gainCurvePower2(hid2, semilog=False, figp=False, linestyle='--',label='rm chirp')
plt.legend()
plt.show()
stp.savefig(filename="compare_energy")

#%%% compare, rotate undulator 90 deg
from xkinterface.bbl import postgen4

fname='/mnt/f/simu_2025/202510_idealMachineTHzFEL/0000_syncDESY/beamEnergy_17MeV/gauss_Q2nC_phiboo_-40deg/02_gen4/exchange_xy_compress_in_xplane/gen4.out.h5'
hid = h5py.File(fname,'r')

fname='/mnt/f/simu_2025/202510_idealMachineTHzFEL/0000_syncDESY/beamEnergy_17MeV/gauss_Q2nC_phiboo_-40deg/02_gen4/nocollim/gen4.out.h5'
hid2 = h5py.File(fname,'r')

fig,ax1 = plt.subplots()
postgen4.plot_gainCurvePeakPower2(hid, semilog=False, figp=False,linestyle='-',label='kx=1,ky=0')

postgen4.plot_gainCurvePeakPower2(hid2,semilog=False, figp=False,linestyle='--',label='kx=0,ky=1')
plt.legend()

stp.savefig('gain_compare')


plt.figure()
postgen4.plot_gainCurvePower2(hid, semilog=False, figp=False, linestyle='-',label='kx=1,ky=0')
postgen4.plot_gainCurvePower2(hid2, semilog=False, figp=False, linestyle='--',label='kx=0,ky=1')
plt.legend()
plt.show()
stp.savefig(filename="compare_energy")


#%% rms size
from xkinterface.bbl import postgen4

fname = 'gen4.out.h5'
hid = h5py.File(fname,'r')


# postgen4.plot_rmsSizeEvo(hid, field=False)
# stp.savefig(filename="rmsSize")

# postgen4.plot_rmsSizeEvo(hid, field=True)
# plt.grid(True)
# stp.savefig(filename="rmsSize_field")


# postgen4.plot_rmsSizeEvo(hid, field=False, fig=False)

#%%% rms size of certain slice
from xkinterface.bbl import postgen4

postgen4.plot_rmsSizeEvo_sliceJ(hid, slicej=30)

#%%% beam center offset
plt.figure()

z = hid['Lattice']['zplot'][()]

x = hid["Beam"]["Global"]["xposition"][:]
y = hid["Beam"]["Global"]["yposition"][:]

plt.plot(z, x*1e3,'--',label='x0')
plt.plot(z,y*1e3, '-', label='y0')

plt.xlabel("s (m)")
plt.ylabel("beam center (mm)")
plt.legend()

#%% power profile evo
from xkinterface.interface import postgen4

# path='/mnt/f/simu_2025/202508_genesis_waveguide_weiweiLi/0915_BCON_gaussian1nC/waveguide_off/gauss_1nC_BCON_phi2-32deg_scOFF'
# os.chdir(path)

fname = "gen4.out.h5"
hid = h5py.File(fname,'r')

postgen4.plot_powerEvolution(hid, norm=False, slicenum=True)
stp.savefig('powerProfielEvo')

postgen4.plot_powerEvolution(hid, norm=True, slicenum=True)
stp.savefig("powerProfile")

#%%% profile line
# pos_sat = 0.8
# jj = int(pos_sat/0.015)

pos_sat = 5.4
jj = int(pos_sat/(110e-3/2))

s = hid['Global']['s'][()]/2.998e8*1e12  #ps
power =  hid['Field']['power'][jj,:]
plt.figure()
plt.plot(s,power/1e6,'-')
plt.xlabel('z (ps)')
plt.ylabel('power (MW)')
plt.title(f's={pos_sat:.2f} m')
stp.savefig(f'power_profile_{pos_sat:.2f}m')


#%%%% compare
pos_sat = 5
half_period=110e-3/2

jj = int(pos_sat/half_period)

s = hid['Global']['s'][()]/2.998e8*1e12  #ps
power =  hid['Field']['power'][jj,:]
plt.figure()
plt.plot(s,power/1e6,'-',label='keep chirp')

s = hid2['Global']['s'][()]/2.998e8*1e12  #ps
power =  hid2['Field']['power'][jj,:]
plt.plot(s,power/1e6,'--',label='rm chirp')

plt.xlabel('z (ps)')
plt.ylabel('power (MW)')
plt.legend()
plt.title(f's={pos_sat:.2f} m')
stp.savefig(f'power_profile_{pos_sat:.2f}m')

#%%%% compare-2
pos_sat = 0.8
half_period=30e-3/2

jj = int(pos_sat/half_period)

s = hid['Global']['s'][()]/2.998e8*1e12  #ps
power =  hid['Field']['power'][jj,:]
plt.figure()
plt.plot(s,power/1e6,'-',label='no alignment')

s = hid2['Global']['s'][()]/2.998e8*1e12  #ps
power =  hid2['Field']['power'][jj,:]
plt.plot(s,power/1e6,'--',label='alignment')

plt.xlim([-1, 20])

plt.xlabel('z (ps)')
plt.ylabel('power (MW)')
plt.legend()
plt.title(f's={pos_sat:.2f} m')
stp.savefig(f'power_profile_{pos_sat:.2f}m')

#%%%
fname = './gen4.out.h5'
hid = h5py.File(fname,'r')

steps_num = hid["Field"]["power"].shape[0]
jjl = np.arange(0,steps_num,5)

for jj in jjl:
    # spos = 0.015*jj*3
    spos = 0.015*jj*1

    s = hid['Global']['s'][()]/2.998e8*1e12  #ps
    power =  hid['Field']['power'][jj,:]
    
    plt.figure()
    plt.plot(s,power/1e6,'-')
    plt.xlabel('z (ps)')
    plt.ylabel('power (MW)')
    plt.title(f's={spos:.2f} m')
    stp.savefig(f'power_profile_{spos:.2f}m')
    plt.close()

#%% beam slice energy evo
from xkinterface.interface import postgen4

postgen4.plot_beamEnergyEvolution(hid, sliceRange=[0,85])
stp.savefig("sliceEnergy")

#%% beam slice bunching evo
from xkinterface.interface import postgen4

postgen4.plot_beamBunchingEvolution(hid, norm=True) #, sliceRange=[5,90])
stp.savefig("sliceBunchingNorm")

postgen4.plot_beamBunchingEvolution(hid, norm=False) #, sliceRange=[5,90])
stp.savefig("sliceBunching")
#%% spectrum
from xkinterface.bbl import postgen4

postgen4.plot_spectrum(hid, step=-1)
stp.savefig("spectrum")

# postgen4.plot_spectrumEvolution(hid, norm=False) #, freqRange=[2.5,3.5])
# stp.savefig("spectrumEvo")

# postgen4.plot_spectrumEvolution(hid, norm=True) #, freqRange=[2.5,3.5])
# stp.savefig("spectrumEvo_norm")

#%% Beam phase space evolution
# from xkinterface.interface import Plot as xkplt
import impzpy.post as pzplot
import h5py

sampleFreq=1

step = [int(0.1/0.15)]
# step = [int(3.0/0.15)]
# step = [int(2.25/0.15)]

jjl = np.arange(0,260,20)
# jjl=[1,2]

# for jj in range(25):
for jj in jjl:
    # fname = f"dump{jj}.par.h5"
    fname = f'gen4.{jj}.par.h5'

    f = h5py.File(fname,'r')

    position_s = jj*0.015 #m
    
    # lamdas = 15e-6
    lamdas = 100e-6
    
    s1 = 0
    s2 = 100
    dlamda = s1
    
    x= np.array([])
    y= np.array([])
    # for j in range(90):
    for j in np.arange(s1,s2,1):

        name = f'slice{j+1:06d}'
        
        x = np.append(x,f[name]['theta'][:]/(2*np.pi)+dlamda)
        y = np.append(y,f[name]['gamma'][:])
        # plt.plot(f[name]['theta'][:]/np.pi*180+dphi,f[name]['gamma'][:],'.r')
        dlamda = dlamda +1
    
    x = x-np.min(x)
    
    plt.figure(figsize=(10*1.2,4*1.2))
    
    plt.suptitle(f"s={position_s:.2f} m")
    plt.subplot(1,2,1)
    # plt.title(f"@{position_s:.2f} m")
    plt.plot(x[::sampleFreq],y[::sampleFreq],'r.', markersize=1)
    plt.xlabel(r'$z/\lambda_s$')
    plt.ylabel(r'$\gamma$')
    plt.grid()

    plt.subplot(1,2,2)
    # plt.plot(x[::sampleFreq],y[::sampleFreq],'r.', markersize=1)
    
    x1=40; x2=45
    tmp1 = x[(x>x1) & (x<x2)]
    tmp2 = y[(x>x1) & (x<x2)]
    # plt.plot(tmp1[::sampleFreq],tmp2[::sampleFreq],'r.', markersize=1)
    pzplot.plot_dotpha(tmp1,tmp2,bins=(x2-x1)*10)
    # plt.ylim([28,36])
    plt.xlabel(r'$z/\lambda_s$')
    plt.ylabel(r'$\gamma$')
    plt.grid()
    

    plt.savefig(f"{position_s:.2f}m_beamPhase_22.png")
    
    if len(jjl) > 1:
        plt.close()


    # ##%% current from dumped beam
    # bins =256*4
    # lamdas = 100e-6   #15e-6
    
    # zz = x*lamdas
    
    
    
    # hist, bin_edges = np.histogram(zz,bins=bins)
    
    # dt = (bin_edges[2]-bin_edges[1])/2.998e8
    # charge = np.sum( hid["Beam"]["current"][0,:] ) *lamdas/2.998e8
    # qi = charge/len(zz)
    # Ii = qi/dt
    
    # z = bin_edges[0:-1]/2.998e8 *1e12  #ps
    
    # plt.figure()    
    # plt.title(f"s={position_s:.2f} m")    
    # plt.plot(z,hist*Ii,'b-')
    # plt.xlabel("z (ps)")
    # plt.ylabel("current (A)")
    # # plt.grid()
    
    # stp.savefig("current")

# ##%% current from dumped beam
# bins =256*4
# lamdas = 100e-6   #15e-6

# zz = x*lamdas



# hist, bin_edges = np.histogram(zz,bins=bins)

# dt = (bin_edges[2]-bin_edges[1])/2.998e8
# charge = np.sum( hid["Beam"]["current"][0,:] ) *lamdas/2.998e8
# qi = charge/len(zz)
# Ii = qi/dt

# z = bin_edges[0:-1]/2.998e8 *1e12  #ps

# plt.figure()    
# plt.title(f"s={position_s:.2f} m")    
# plt.plot(z,hist*Ii,'b-')
# plt.xlabel("z (ps)")
# plt.ylabel("current (A)")
# # plt.grid()

# stp.savefig("current")


#%%% single phase space
import impzpy.post as pzplot

fname = f'inibeam.par.h5'

f = h5py.File(fname,'r')

# lamdas = 15e-6
lamdas = 60e-6
sampleFreq = 1

s1 = 0
s2 = 180
dlamda = s1

x= np.array([])
y= np.array([])
# for j in range(90):
for j in np.arange(s1,s2,1):

    name = f'slice{j+1:06d}'
    
    x = np.append(x,f[name]['theta'][:]/(2*np.pi)+dlamda)
    y = np.append(y,f[name]['gamma'][:])
    # plt.plot(f[name]['theta'][:]/np.pi*180+dphi,f[name]['gamma'][:],'.r')
    dlamda = dlamda +1

x = x-np.min(x)

plt.figure(figsize=(10*1.2,4*1.2))

plt.subplot(1,2,1)
# plt.title(f"@{position_s:.2f} m")
plt.plot(x[::sampleFreq],y[::sampleFreq],'r.', markersize=1)
plt.xlabel(r'$z/\lambda_s$')
plt.ylabel(r'$\gamma$')
plt.grid()

# plt.subplot(1,2,2)
# plt.plot(x[::sampleFreq],y[::sampleFreq],'r.', markersize=1)
plt.figure()
x1=0; x2=10
tmp1 = x[(x>x1) & (x<x2)]
tmp2 = y[(x>x1) & (x<x2)]
# plt.plot(tmp1[::sampleFreq],tmp2[::sampleFreq],'r.', markersize=1)
pzplot.plot_dotpha(tmp1,tmp2,bins=(x2-x1)*10)
# plt.ylim([28,36])
plt.xlabel(r'$z/\lambda_s$')
plt.ylabel(r'$\gamma$')
plt.grid()

#%% field plot
import os
import matplotlib.pyplot as plt
import h5py 

path='/mnt/f/simu_2025/202508_genesis_waveguide_weiweiLi/0930_pitz_paras_PPS_Pre/Gen4_Waveguide/gen4_waveguide'
os.chdir(path)

fname = "gen4.out.h5"
hid = h5py.File(fname,'r')

def getWF(filename,slice=1):
    hfl = h5py.File(filename,'r')
    slc = 'slice%6.6d' % slice
    ng = hfl['gridpoints'][()][0]
    dg = hfl['gridsize'][()][0]
    fre = hfl[slc]['field-real'][()]
    fim = hfl[slc]['field-imag'][()]
    inten = np.reshape(fre*fre+fim*fim, (ng,ng))
    return np.sqrt(inten),dg*(ng-1)*0.5*1e3, hfl


# stepjl = np.arange(0,240,20)
slil = np.arange(0,150,1)

stepjl=[226]

for stepj in stepjl:        
    maxPower_sli_id = np.argmax(hid['Field']['power'][stepj,:])
    
    intentol = np.zeros((301,301))
    for jj in slil:
    # for jj in [maxPower_sli_id]:
        fname = f"gen4.{stepj}.fld.h5"
        
        inten, dg, hfl = getWF(fname,slice=jj+1)
        
        intentol += inten
    
        # plt.figure()
        # plt.title(f"slice {jj}")
        
        # # inten /= np.max( inten )
        # masked_profile = np.where(inten < 0.001, np.nan, inten)
        
        # plt.imshow(masked_profile,extent=(-dg,dg,-dg,dg),vmin=0.001)
        # plt.colorbar()
    plt.figure()
    plt.title(f"slice {jj}")
        
    # inten /= np.max( inten )
    # masked_profile = np.where(inten < 0.001, np.nan, intentol)
        
    plt.imshow(intentol,extent=(-dg,dg,-dg,dg),vmin=0.001)
    plt.colorbar()

#%%% field plot
def getWF(filename,slice=1):
    hfl = h5py.File(filename,'r')
    slc = 'slice%6.6d' % slice
    ng = hfl['gridpoints'][()][0]
    dg = hfl['gridsize'][()][0]
    fre = hfl[slc]['field-real'][()]
    fim = hfl[slc]['field-imag'][()]
    inten = np.reshape(fre*fre+fim*fim, (ng,ng))
    return np.sqrt(inten),dg*(ng-1)*0.5*1e3, hfl


path='/mnt/f/simu_2025/202508_genesis_waveguide_weiweiLi/0930_pitz_paras_PPS_Pre/Gen4_Waveguide/gen4_waveguide'
os.chdir(path)

fname = "gen4.out.h5"
hid = h5py.File(fname,'r')

slil = np.arange(0,150,1)
intentol = np.zeros((89,89))
for jj in slil:
# for jj in [maxPower_sli_id]:
    fname = f"gen4.226.fld.h5"
    
    inten, dg, hfl = getWF(fname,slice=jj+1)
    
    intentol += inten

    # plt.figure()
    # plt.title(f"slice {jj}")
    
    # # inten /= np.max( inten )
    # masked_profile = np.where(inten < 0.001, np.nan, inten)
    
    # plt.imshow(masked_profile,extent=(-dg,dg,-dg,dg),vmin=0.001)
    # plt.colorbar()
plt.figure()
plt.title(f"slice {jj}")
    
# inten /= np.max( inten )
# masked_profile = np.where(inten < 0.001, np.nan, intentol)
    
plt.imshow(intentol,extent=(-dg,dg,-dg,dg),vmin=0.001,cmap='viridis')
plt.colorbar()


#%% beam phase space
import impzpy.post as pzplot
# path = stp.winpath(r'F:\simu_2025\202509_XiangkunSimu\beam_112A_2.0nC_shotnoise_v4_quiet2')
# os.chdir(path)

fname = r'368A.2809.002.1000.h5'
fname = r'368A.2809.002.1000_5.4m.h5'
fname=r'dump.par.h5'
fname=r'gen4_sliced_one4one.h5'


f = h5py.File(fname,'r')

# lamdas = 15e-6
lamdas = 100e-6
sampleFreq = 10

s1 = 0
s2 = 99
dlamda = s1

x= np.array([])
y= np.array([])
# for j in range(90):
for j in np.arange(s1,s2,1):

    name = f'slice{j+1:06d}'
    
    x = np.append(x,f[name]['theta'][:]/(2*np.pi)+dlamda)
    y = np.append(y,f[name]['gamma'][:])
    # plt.plot(f[name]['theta'][:]/np.pi*180+dphi,f[name]['gamma'][:],'.r')
    dlamda = dlamda +1

x = x-np.min(x)

plt.figure(figsize=(10*1.2,4*1.2))

# plt.suptitle(f"s={position_s:.2f} m")
plt.subplot(1,2,1)
# plt.title(f"@{position_s:.2f} m")ls
plt.plot(x[::sampleFreq],y[::sampleFreq],'r.', markersize=1)
plt.xlabel(r'$z/\lambda_s$')
plt.ylabel(r'$\gamma$')
plt.grid()

plt.subplot(1,2,2)
# plt.plot(x[::sampleFreq],y[::sampleFreq],'r.', markersize=1)

# x1=40; x2=45
x1=0; x2=80
tmp1 = x[(x>x1) & (x<x2)]
tmp2 = y[(x>x1) & (x<x2)]
# plt.plot(tmp1[::sampleFreq],tmp2[::sampleFreq],'r.', markersize=1)
pzplot.plot_dotpha(tmp1,tmp2,bins=(x2-x1)*10)
# plt.ylim([28,36])
plt.xlabel(r'$z/\lambda_s$')
plt.ylabel(r'$\gamma$')
plt.grid()

stp.savefig('beamphase')

#%%% x-y phase space
import impzpy.post as pzplot
# path = stp.winpath(r'F:\simu_2025\202509_XiangkunSimu\beam_112A_2.0nC_shotnoise_v4_quiet2')
# os.chdir(path)
fname=r'gen4_sliced_one4one.h5'

f = h5py.File(fname,'r')

# lamdas = 15e-6
lamdas = 100e-6
sampleFreq = 1

s1 = 0
s2 = 99
dlamda = s1

x= np.array([])
y= np.array([])
# for j in range(90):
for j in np.arange(s1,s2,1):

    name = f'slice{j+1:06d}'
    
    x = np.append(x,f[name]['x'][:])
    y = np.append(y,f[name]['y'][:])
    # plt.plot(f[name]['theta'][:]/np.pi*180+dphi,f[name]['gamma'][:],'.r')

plt.figure(figsize=(10*1.2,4*1.2))

# plt.title(f"@{position_s:.2f} m")ls
plt.plot(x[::sampleFreq],y[::sampleFreq],'r.', markersize=1)
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.grid()

#%%% x-xp,,, 6D phase
import impzpy.post as pzplot
path = stp.winpath(r'F:\simu_2025\initialSamplingMethod\plannerUndu\04_xyFocus_smoothGaussProfile\test_timeON_xyFocus')

path='/mnt/f/simu_2025/initialSamplingMethod/plannerUndu/04_xyFocus_smoothGaussProfile/00_genDis_impz'
os.chdir(path)



# fname = r'368A.2809.002.1000.h5'
fname = r'gen4_sliced_one4one.h5'
# fname = r'368A.2809.002.1000_5.4m.h5'
# fname=r'dump1.par.h5'

f = h5py.File(fname,'r')

xdat = f['theta'][:]
ydat = f['x'][:]

plt.figure()
pzplot.plot_dotpha(xdat, ydat,'.')















