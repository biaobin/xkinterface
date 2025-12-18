#%% initial setting

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from impzpy import post as impz_post
from scipy.interpolate import interp1d
import h5py
import xkinterface.interface.PostGenesis13 as PostGenesis
import sys
from startup import *
    
# workdir='/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_2nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_-24.0deg_phi1_5deg_sigt4ps/02_gen4'
# workdir='/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/Q_2nC_1Len_lscOFF'
# workdir="/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_4nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_-20.0deg_phi1_0deg"
# workdir="/lustre/fs22/group/pitz/biaobin/2025/flattop_bcON_4nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_-16.0deg_phi1_0deg"
workdir="/mnt/f/simu_2025/test_100MeV"
workdir = "/mnt/f/simu_2025/test_17MeV_one4one_scOFF"
workdir = "/home/ubuntu2404/github/Genesis-1.3-Version4/examples/Example1-SteadyState"

workdir = winpath(r"F:\simu_2025\initialSamplingMethod\plannerUndu\04_match2twiss_timeON_SmoothGaussProfile\test_timeOFF")

# fname = 'Example1.out.h5'
fname = "gen4.out.h5"
os.chdir(workdir)    
    

#%% impz section
from impzpy import post as impz_post
path = "./01_impz/"

#%%% rms & beta function plot 
twi = impz_post.twiss(path=path)

twi.plot_twiss()
twi.plot_rms(enxy="ON") 

#%%% phase space: x-xp, x-y, z-dgam
pha = impz_post.phase(fname="fort.1000",path=path)

pha.plot_dot_xxp_xy_zdE(sample=10)
# pha.plot_heat_xxp_xy_zdE()

# plt.figure()
# plt.hist2d(pha.z,pha.dgam,bins=64)

#%%% long. phase space evolution
from impzpy import post as impz_post
path = "./01_impz/"

deg = -14
workdir=f"/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_4nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_{deg}.0deg_phi1_0deg"
workdir=f"/lustre/fs22/group/pitz/biaobin/2025/flattop_bcON_4nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_{deg}.0deg_phi1_0deg"
workdir=f"/lustre/fs22/group/pitz/biaobin/2025/flattop_bcON_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{deg}.0deg"

os.chdir(workdir)  

plt.figure(figsize=(14,6))
plt.suptitle(f"phi2={deg} deg")
pha = impz_post.phase(fname="fort.1001",path=path)
pha.plot_zdgam_current(dot=True, linelabel="before BC",profile=True )

pha = impz_post.phase(fname="fort.1002",path=path)
pha.plot_zdgam_current(dot=True, linelabel="after BC",profile=True)  

# pha = impz_post.phase(fname="fort.1003",path=path)
# pha.plot_zdgam_current(dot=True)    

#%%%% longi phase for flattop-2nC case
from impzpy import post as impz_post
path = "./01_impz/"

deg = -16
workdir=f"/lustre/fs22/group/pitz/biaobin/2025/flattop_bcON_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{deg}.0deg"

os.chdir(workdir)  

plt.figure(figsize=(7,6))
pha = impz_post.phase(fname="fort.1000",path=path)
pha.plot_zdgam(dot=True, linelabel="boos. exit",profile=True)  
plt.figure(figsize=(7,6))
pha = impz_post.phase(fname="fort.1001",path=path)
pha.plot_zdgam(dot=True, linelabel="before BC",profile=True)  
plt.figure(figsize=(7,6))
pha = impz_post.phase(fname="fort.1002",path=path)
pha.plot_zdgam(dot=True, linelabel="after BC",profile=True)  
plt.figure(figsize=(7,6))
pha = impz_post.phase(fname="fort.1003",path=path)
pha.plot_zdgam(dot=True, linelabel="undu. entr.",profile=True)  

# # BC OFF
# #================
# deg = -10
# workdir=f"/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{deg}.0deg"

# os.chdir(workdir)  

# plt.figure(figsize=(7,6))
# pha = impz_post.phase(fname="fort.1000",path=path)
# pha.plot_zdgam(dot=True, linelabel="boos. exit",profile=True)  

# plt.figure(figsize=(7,6))
# pha = impz_post.phase(fname="fort.1003",path=path)
# pha.plot_zdgam(dot=True, linelabel="undu. entr.",profile=True)  

#%%%% longi phase for 1nC, gauss
# for 2025-collab-meeting
from impzpy import post as impz_post
from xkinterface.interface import Plot as xkplt

workdir='/lustre/fs22/group/pitz/biaobin/2025/gauss_bcOFF_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/01_impz'
os.chdir(workdir)  

# pha = impz_post.phase(fname="fort.1000")
# plt.figure(figsize=(8,6.2))
# pha.plot_zdgam(dot=True, profile=True, bins=64,linelabel=r"Booster exit")

# pha = impz_post.phase(fname="fort.1001")
# plt.figure(figsize=(8,6.2))
# pha.plot_zdgam(dot=True, profile=True, bins=64,linelabel=r"BC entr.")

# pha = impz_post.phase(fname="fort.1002")
# plt.figure(figsize=(8,6.2))
# pha.plot_zdgam(dot=True, profile=True, bins=64,linelabel=r"BC exit")

pha = impz_post.phase(fname="fort.1003")
plt.figure(figsize=(8,6.2))
# pha.plot_zdgam_current(dot=True, profile=True, bins=64,linelabel=r"undu. entr.")
pha.plot_zdgam(dot=True, profile=True, bins=100)
# pha = impz_post.phase(fname="fort.1001")
# plt.figure(figsize=(8,6.2))
# pha.plot_zdgam(dot=True, profile=True, bins=64,linelabel=r"BC entr.")

# pha = impz_post.phase(fname="fort.1002")
# plt.figure(figsize=(8,6.2))4,linelabel=r"undu. entr.")

# pha.plot_zdgam_current(dot=True, linelabel="before BC")
# xkplt.plot2d(pha.zt,pha.dgam,xlabel="z (ps)",ylabel=r"$\Delta \gamma$")

# plt.figure()
# ax,ay, sgnl = xkplt.hist2d_contourf(pha.zt, pha.dgam)
# sgnl = abs(sgnl)
# sgnl = sgnl/np.max(sgnl)
# v = np.linspace(0.001, 1, 99)
# plt.contourf(ax,ay,sgnl,v)


#%%%% test
#%%%% longi phase for flattop-2nC case
from impzpy import post as impz_post
path = "./01_impz/"

deg = -16
workdir=f"/lustre/fs22/group/pitz/biaobin/2025/flattop_bcON_2nC/cores_016_impt_
# ax,ay, sgnl = xkplt.hist2d_contourf(pha.zt, pha.dgam)64640064_impz_64640064_Np_5.24e+06_phi2_{deg}.0deg"

os.chdir(workdir)  

plt.figure(figsize=(7,6))
pha = impz_post.phase(fname="fort.1000",path=path)
pha.plot_zdgam(dot=True, linelabel="boos. exit",profile=True)  
plt.figure(figsize=(7,6))
pha = impz_post.phase(fname="fort.1001",path=path)
pha.plot_zdgam(dot=True, linelabel="before BC",profile=True)  
plt.figure(figsize=(7,6))
pha = impz_post.phase(fname="fort.1002",path=path)
pha.plot_zdgam(dot=True, linelabel="after BC",profile=True)  
plt.figure(figsize=(7,6))
pha = impz_post.phase(fname="fort.1003",path=path)
pha.plot_zdgam(dot=True, linelabel="undu. entr.",profile=True)  




#%% genesis section

#%%% gain curve
import xkinterface.interface.PostGenesis13 as PostGenesis

workdir="/lustre/fs22/group/pitz/biaobin/2025/flattop_bcON_4nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_-16.0deg_phi1_0deg"
# workdir="/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_4nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_-20.0deg_phi1_0deg"
workdir= f"/lustre/fs22/group/pitz/biaobin/2025/gauss_bcOFF_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg"
os.chdir(workdir) 

fname = "./02_gen4/gen4.seed5.nper10.fitdeg3.out.h5"
fname = "./02_gen4/gen4.seed5.nper5.fitdeg1.out.h5"

workdir=f"/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4/slsc_off_lscON"
os.chdir(workdir) 
fname = "gen4.out.h5"

pg1 = PostGenesis(fname, debug=1, version=4, fig_ext = '.png')

# pg1.plot_tpower(at=2.6)
# pg1.plot_spectrum(at=2.6)
# pg1.plot_current()

# h5f = h5py.File(fname,'r')

# plt.figure()
# plt.plot(pg1.zplot,pg1.zenergy*1e6,'-')    

#%%% check the slippage effects
plt.figure(figsize=(8,10))
plt.subplot(2,1,1)
pg1.plot_current(x="SLICE",fig=False)
plt.subplot(2,1,2)
pg1.plot_tpower(x="SLICE",fig=False, at=1.5)


#%%% compare gain curve 
#======================


#%%%% cut at two sides
workdir='/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-6.0deg/02_gen4'
os.chdir(workdir)
fname = "./gen4.seed5.nper5.fitdeg1.out.h5"
pg0 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')


workdir='/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/'
workdir='/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-6.0deg/02_gen4/cut_at_two_side_sllscOFF'
os.chdir(workdir)

fname = "./gen4.seed5.nper5.fitdeg1.out.h5"
pg1 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

# pg1.plot_tpower(at=2.6)
# pg1.plot_spectrum(at=2.6)
# pg1.plot_current()


workdir='/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-30.0deg/02_gen4'
workdir='/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_2nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_-28.0deg_phi1_5deg_sigt4ps/02_gen4'
workdir='/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/Q_1nC_1Len_lscOFF'
workdir='/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-6.0deg/02_gen4/cut_at_two_side'
os.chdir(workdir)

fname = "./gen4.seed5.nper5.fitdeg1.out.h5"
pg2 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

plt.figure(figsize=(8,7))
pg0.plot_energy(fig=False, line="-", label="no cut, S-L-LSC ON")
pg2.plot_energy(fig=False, line="-.", label="cut, S-L-LSC ON")
pg1.plot_energy(fig=False, line="--", label="cut, S-L-LSC OFF")
plt.grid()

plt.figure(figsize=(8,7))
pg0.plot_tpower(fig=False, linestyle="-", linelabel="no cut, S-L-LSC ON", at=3.0)
pg2.plot_tpower(fig=False, linestyle="-.", linelabel="cut, S-L-LSC ON",   at=3.0)
pg1.plot_tpower(fig=False, linestyle="--", linelabel="cut, S-L-LSC OFF",  at=2.5)
plt.grid()


# pg2.plot_tpower(at=2.6)
# pg2.plot_spectrum(at=2.6)
# pg2.plot_current()
   
# plt.figure()
# plt.plot(pg1.zplot,pg1.zenergy*1e6,'-',label="phi1_5_phi2_-24")    
# plt.plot(pg2.zplot,pg2.zenergy*1e6,'--',label="phi1_5_phi2_-28")    
# plt.yscale("log")
# plt.xlabel("s (m)")
# plt.ylabel("THz energy (uJ)")
# plt.grid()
# plt.legend()

#%%%% add peak at tail
workdir='/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-6.0deg/02_gen4'
os.chdir(workdir)
fname = "./gen4.seed5.nper5.fitdeg1.out.h5"
pg0 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')


workdir='/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-6.0deg/02_gen4/addPeak'
os.chdir(workdir)
fname = "./gen4.seed5.nper10.fitdeg3.out.h5"
pg1 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')


workdir='/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-6.0deg/02_gen4/addPeak_sllscOFF'
os.chdir(workdir)

fname = "./gen4.seed5.nper10.fitdeg3.out.h5"
pg2 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

plt.figure(figsize=(8,7))
pg0.plot_energy(fig=False, line="-", label="no peak, S-L-LSC ON")
pg1.plot_energy(fig=False, line="-.", label="add peak, S-L-LSC ON")
pg2.plot_energy(fig=False, line="--", label="add peak, S-L-LSC OFF")
plt.grid(True)

plt.figure(figsize=(8,7))
pg0.plot_tpower(fig=False, linestyle="-", linelabel="no peak, S-L-LSC ON", at=3.0)
pg1.plot_tpower(fig=False, linestyle="-.", linelabel="add peak, S-L-LSC ON",   at=3.0)
pg2.plot_tpower(fig=False, linestyle="--", linelabel="add peak, S-L-LSC OFF",  at=3.0)
plt.grid(True)


#%%%% larger beam size
workdir="/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4"
os.chdir(workdir)
fname = "./gen4.seed5.nper10.fitdeg3.out.h5"
pg0 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

workdir="/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/test_1THz"
os.chdir(workdir)
fname = "./gen4.seed5.nper10.fitdeg3.out.h5"
pg1 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')



plt.figure(figsize=(8,7))
pg0.plot_energy(fig=False, line="-", label="sigx=0.5mm, S-L-LSC ON")
pg1.plot_energy(fig=False, line="-.", label="sigx=1mm, S-L-LSC ON")
plt.grid(True)

plt.figure(figsize=(8,7))
pg0.plot_tpower(fig=False, linestyle="-", linelabel="sigx=0.5mm, S-L-LSC ON")
pg1.plot_tpower(fig=False, linestyle="-.", linelabel="sigx=1mm, S-L-LSC ON")
plt.grid(True)

#%%%% hollow beam
label_ho="hollow, rb=0.3 mm"

workdir="/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4"
os.chdir(workdir)
fname = "./gen4.seed5.nper10.fitdeg3.out.h5"
pg0 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

workdir="/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/test_hollowBeam"
os.chdir(workdir)
fname = "./gen4.seed5.nper10.fitdeg3.out.h5"
pg1 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

plt.figure(figsize=(8,7))
pg0.plot_energy(fig=False, line="-", label="no hollow, S-L-LSC ON")
pg1.plot_energy(fig=False, line="-.", label=f"{label_ho}, S-L-LSC ON")
plt.grid(True)

plt.figure(figsize=(8,7))
pg0.plot_tpower(fig=False, linestyle="-", linelabel="no hollow, S-L-LSC ON")
pg1.plot_tpower(fig=False, linestyle="-.", linelabel=f"{label_ho}, S-L-LSC ON")
plt.grid(True)

plt.figure(figsize=(8,7))
pg0.plot_spectrum(fig=False, linestyle="-", linelabel="no hollow, S-L-LSC ON")
pg1.plot_spectrum(fig=False, linestyle="-.", linelabel=f"{label_ho}, S-L-LSC ON")
plt.grid(True)

#%%%% 50 MeV beam
label_ho="50 MeV"

# workdir="/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4"
workdir="/mnt/f/simu_2025/test_50MeV/test_nzlsc_1"
os.chdir(workdir)
fname = "./gen4.seed5.nper10.fitdeg3.out.h5"
pg0 = PostGenesis(fname, debug=1, version=4, fig_ext = '.png')

# workdir="/afs/ifh.de/user/b/biaobin/sim2/2025/test_50MeV"
workdir="/mnt/f/simu_2025/test_50MeV"
os.chdir(workdir)
fname = "./gen4.seed5.nper10.fitdeg3.out.h5"
pg1 = PostGenesis(fname, debug=1, version=4, fig_ext = '.png')

plt.figure(figsize=(8,7))
pg0.plot_energy(fig=False, line="-", label="no hollow, S-L-LSC ON")
pg1.plot_energy(fig=False, line="-.", label=f"{label_ho}, S-L-LSC ON")
plt.grid(True)

plt.figure(figsize=(8,7))
pg0.plot_tpower(fig=False, linestyle="-", linelabel="no hollow, S-L-LSC ON")
pg1.plot_tpower(fig=False, linestyle="-.", linelabel=f"{label_ho}, S-L-LSC ON")
plt.grid(True)

plt.figure(figsize=(8,7))
pg0.plot_spectrum(fig=False, linestyle="-", linelabel="no hollow, S-L-LSC ON")
pg1.plot_spectrum(fig=False, linestyle="-.", linelabel=f"{label_ho}, S-L-LSC ON")
plt.grid(True)

#%%%% 25-04-10, for report, gauss 1nC, slsc, llsc
workdir=f"/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4"
os.chdir(workdir)  
fname = "gen4.out.h5"
pg1 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

workdir=f"/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4/slsc_off_lscON"
os.chdir(workdir)
fname = "gen4.out.h5"
pg2 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

workdir=f"/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4/slsc_off"
os.chdir(workdir)
fname = "gen4.out.h5"
pg3 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

workdir=f"/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/20241202-NonCompressor-ScanPhase/cores_016_impt_32320064_impz_32320064_Np_1.00e+07_phi2_-36.0deg/02_gen4"
os.chdir(workdir)
fname = "gen4.out.h5"
pg4 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')


plt.figure(figsize=(7,6))
plt.plot(pg1.zplot,pg1.zenergy*1e6,'-',label="BC ON, short-long-range LSC ON")      
plt.plot(pg3.zplot,pg3.zenergy*1e6,'--',label="BC ON, short-long-range LSC OFF")  
plt.plot(pg2.zplot,pg2.zenergy*1e6,'-.',label="BC ON, short-range LSC OFF")  
plt.plot(pg4.zplot,pg4.zenergy*1e6,'-',label="BC OFF, short-long-range LSC ON")   
plt.yscale("log")
plt.xlabel("s (m)")
plt.ylabel("THz energy (uJ)")
plt.grid()
plt.legend()

#%%%% 25-04-10, for report, t-power profile
import xkinterface.interface.PostGenesis13 as PostGenesis

#======================
plt.figure(figsize=(7,6))
pg1.plot_tpower(fig=False, at=2.0, linelabel="BC ON, S-L-LSC ON",linestyle="-")
pg3.plot_tpower(fig=False, at=1.8, linelabel="BC ON, S-L-LSC OFF",linestyle="--")
pg2.plot_tpower(fig=False, at=1.8, linelabel="BC ON, S-LSC OFF",linestyle="-.")
pg4.plot_tpower(fig=False, at=2.4, linelabel="BC OFF, S-L-LSC ON",linestyle="-", xrange=[0,30])
plt.grid()


#%%%% 25-04-21, for report, spectrum
import xkinterface.interface.PostGenesis13 as PostGenesis

plt.figure(figsize=(7,6))
pg1.plot_spectrum(fig=False, at=2.0, linelabel="BC ON, S-L-LSC ON",linestyle="-", norm=True)
pg3.plot_spectrum(fig=False, at=1.8, linelabel="BC ON, S-L-LSC OFF",linestyle="--",norm=True)
pg2.plot_spectrum(fig=False, at=1.8, linelabel="BC ON, S-LSC OFF",linestyle="-.", norm=True)
pg4.plot_spectrum(fig=False, at=2.4, linelabel="BC OFF, S-L-LSC ON",linestyle="-", norm=True,loc="upper right")
plt.grid()

#%%%% for Xuduo's case
workdir=f"/lustre/fs22/group/pitz/duoxup/THzSuperRad/genesis/cluster00000001/f=3-Q=0.2-sigz=150.00-sigt=500.00/000copy"
os.chdir(workdir)  
fname = "g4.000.out.h5"
pg1 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

workdir=f"/lustre/fs22/group/pitz/duoxup/THzSuperRad/genesis/cluster00000001/f=3-Q=0.2-sigz=150.00-sigt=500.00"

os.chdir(workdir)


fname = "g4.000.out.h5"
pg2 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

plt.figure(figsize=(7,6))
plt.plot(pg1.zplot,pg1.zenergy*1e6,'--',label="LSC OFF")      
 
plt.plot(pg2.zplot,pg2.zenergy*1e6,'-',label="LSC ON")  
 
plt.yscale("log")
plt.xlabel("s (m)")
plt.ylabel("THz energy (uJ)")
plt.grid()
plt.legend()

plt.figure(figsize=(7,6))
pg1.plot_tpower(fig=False, at=pg1.Lsat, linelabel="LSC OFF",linestyle="--",fig_ext=None)
pg2.plot_tpower(fig=False, at=pg2.Lsat, linelabel="LSC ON",linestyle="-",fig_ext=None)
plt.grid()

# plt.figure(figsize=(7,6))
# pg1.plot_tpower(fig=False, at=-1, linelabel="S-L-LSC OFF",linestyle="-",fig_ext=None)
# pg2.plot_tpower(fig=False, at=-1, linelabel="S-L-LSC ON",linestyle="-.",fig_ext=None)
# plt.grid()

#%%%% flattop-2nC case, for report
workdir=f"/lustre/fs22/group/pitz/biaobin/2025/flattop_bcON_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-16.0deg/02_gen4"
os.chdir(workdir)  
fname = "gen4.seed5.nper10.fitdeg3.out.h5"
pg1 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

workdir=f"/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-10.0deg/02_gen4"
os.chdir(workdir)
fname = "gen4.seed5.nper5.fitdeg1.out.h5"
pg2 = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')


plt.figure(figsize=(7,6))
plt.plot(pg2.zplot,pg2.zenergy*1e6,'-',label="BC OFF")      
plt.plot(pg1.zplot,pg1.zenergy*1e6,'--',label="BC ON")  
 
plt.yscale("log")
plt.xlabel("s (m)")
plt.ylabel("THz energy (uJ)")
plt.grid()
plt.legend()

#%%% beam size evolution
# workdir="/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/test_1THz"
# os.chdir(workdir)
fname = "./gen4.seed5.nper10.fitdeg3.out.h5"
pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

version=4

# Plot electron beam size along the undulator
z = pg.zplot
current = pg.current[1:]

if version < 4: # version 2 or 3
    xx = pg.get_fielddata('xrms')
    yy = pg.get_fielddata('yrms')
    xrms = (xx[:,1:] @ current)/np.sum(current); #print(xrms)
    yrms = (pg.get_fielddata('yrms')[:,1:] @ current)/np.sum(current)
    
    #xrms = np.mean(xx[:,1:], axis = 1); #print(xrms)
    #yrms = np.mean(yy[:,1:], axis = 1); #print(xrms)

else: # version 4
    nn = np.sum(current>0)
    xrms = (pg.get_data('Beam', 'xsize')[:,:nn] @ current[:nn])/np.sum(current)
    yrms = (pg.get_data('Beam', 'ysize')[:,:nn] @ current[:nn])/np.sum(current)

fig, ax = plt.subplots(figsize = (5, 3))

ax.plot(z, xrms*1e3, 'r-')
ax.plot(z, yrms*1e3, 'b-')
ax.plot(z, np.sqrt(xrms*yrms)*1e3, 'g-')

ax.grid()

ax.set_xlabel(r'$z$ (m)')
ax.set_ylabel(r'RMS size (mm)')
ax.legend([r'$x$', r'$y$'])


#%%% animation of FEL process
import xkinterface.interface.PostGenesis13 as PostGenesis

# workdir = f"/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4"  #"/Q_1nC_1Len_lscOFF/"
# # os.chdir(workdir)
# workdir='/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/Q_1nC_2Len_lscON/'

# # workdir="/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_4nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_-20.0deg_phi1_0deg"

# # workdir="/lustre/fs22/group/pitz/biaobin/2025/flattop_bcON_4nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_-16.0deg_phi1_0deg"
# os.chdir(workdir)  
# # # fname = "./02_gen4/gen4.seed5.nper5.fitdeg1.out.h5"
# fname = "gen4.seed5.nper10.fitdeg3.out.h5"

# workdir="/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-6.0deg/02_gen4/cut_at_two_side/.."
# os.chdir(workdir)
# fname="gen4.seed5.nper5.fitdeg1.out.h5"

# workdir="/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-6.0deg/02_gen4/addPeak"

# workdir="/afs/ifh.de/user/b/biaobin/sim2/2025/test_50MeV"
# os.chdir(workdir)
# fname="gen4.seed5.nper10.fitdeg3.out.h5"

# workdir=f"/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4"
# os.chdir(workdir) 
# fname = "gen4.out.h5"

# tmp1="/afs/ifh.de/group/pitz/data/lixiangk/sim3/2023/SASEScan2/beam_200A_2.0nC_shotnoise_p/"
# os.chdir("/lustre/fs22/group/pitz/biaobin/xk_idealGaussFlattop")
# fname = tmp1+"pithz.1.out.h5"
# pg = PostGenesis(fname, debug=1, version=3, fig_ext = '.png')

fname = "gen4.out.h5"
# fname = "Example1.out.h5"
pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

# fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(nrows = 5, figsize = (5, 14))
fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(nrows = 5, figsize = (6, 15))
ax12 = ax1.twinx()
ax32 = ax3.twinx()

# plt.show()

name1 = 'intensity-farfield'
name2 = 'phase-farfield'

ppower = pg.get_fielddata('power')
p_mid = pg.get_fielddata(name1)
phi_mid = pg.get_fielddata(name2)
gamma = pg.get_beamdata('energy')
bunching = pg.get_beamdata('bunching')

xlim1 = 0
xlim2 = 100

step = 5
steps =np.arange(0, ppower.shape[0], step)
# steps = [40]
for i in steps:
    
    # ax1.plot(p_mid[i]/p_mid[i].max(), 'r-')
    # ax12.plot(ppower[i]/ppower[i].max(), 'b-')
    # ax1.plot(p_mid[i]/p_mid[i].max(), 'r-',label='intensity-farfield')
    # ax1.plot(ppower[i]/ppower[i].max(), 'b-',label='power')
    ax1.plot(ppower[i]/1e6, 'b-',label='power')
    
    # title = str.format('z = %.3f m,%0d undu.' % (pg.zplot[i], (pg.zplot[i]-0.105)/0.03))
    title = str.format('z = %.2f m' % (pg.zplot[i]))

    ax1.set_title(title)
    ax1.set_xlabel(r'# of slice')
    ax1.set_ylabel(r'On-axis power (MW)')
    #ax12.set_ylabel(r'Power (au)')
    # ax1.set_ylim(0, 2e7)
    ax1.grid()
    
    ax1.set_xlim(0, xlim2)
    # ax1.set_xlim(0,80)
    # ax1.set_xlim(5,60)
    # ax1.set_xlim(10,50)

    
    # ax1.set_ylim(0., 1.1)
    #ax1.set_yscale('log')
    ax1.legend()
    # ax12.set_ylim(0.0, 1.1)
    #ax12.set_yscale('log')
    
    ax2.plot(gamma[0], 'r-',label='ini.')
    ax2.plot(gamma[i], 'b-',label='step i')
    ax2.set_xlabel('# of slice')
    ax2.set_ylabel(r'$\gamma$')
    # ax2.grid()
    # ax2.set_ylim(32, 35)
    
    ax2.set_xlim(0, xlim2)
    # ax2.set_xlim(0, 125)
    # ax2.set_xlim(10,50)
    # ax2.set_xlim(0,30)
    # ax2.set_xlim(20,100)
    
    ax2.legend()
    
    DE = (gamma[i]-gamma[0]) *pg.current*100e-6/2.998e8/1.6e-19 *0.511e6 *1.6e-19 *1e6 
    DE -= pg.cum_dE_lsc[i] *pg.current*100e-6/2.998e8/1.6e-19 #uJ/slice 
    
    # DE = (gamma[i]-gamma[0])*0.511e3 #keV
    #ax3.plot(DE/np.max(np.abs(DE)), 'r-')
    ax3.plot(DE, 'r-')
    
    ax3.set_xlim(0, xlim2)
    # ax3.set_xlim(10,50)
    # ax3.set_xlim(0,30)
    # ax3.set_xlim(20,100)
    
    # ax3.set_ylim(-2, 2)
    ax3.set_xlabel('# of slice')
    ax3.set_ylabel(r'Energy loss (uJ/slice)')
    ax3.grid(True)
    
    ax32.plot(pg.current, 'k-',label='current prof.')
    #ax32.set_ylabel(r'Current (A)')
    ax32.legend()
    
    ax4.plot(bunching[i], 'r-')
    ax4.set_xlabel('# of slice')
    ax4.set_ylabel('Bunching')
    
    ax4.set_ylim(1e-5, 1)
    ax4.set_yscale('log')
    
    ax4.set_xlim(0, xlim2)
    # ax4.set_xlim(10,50)
    # ax4.set_xlim(0,30)
    # ax4.set_xlim(20,100)
    
    ax4.grid(True)
    
    w, s = pg.get_spectrum(pg.zplot[i])
    #ax5.plot(w*1e6, s/np.max(s), '-')
    ax5.plot(w*1e6, s, '-')
    
    
    # w, s = calcSpectrum(p_mid[i], phi_mid[i])
    # ax5.plot(w*1e6, s/np.max(s), '-')
    
    # ax5.set_xlim(80, 120)
    #ax5.set_ylim(0, 1.1)
    # ax5.set_ylim(0, 2e13)
    
    ax5.set_xlabel(r'Wavelength ($\mu$m)')
    ax5.set_ylabel(r'Intensity (arb. unit)')
    
    # display.display(plt.gcf())
    # display.clear_output(wait=True)

    fig.tight_layout()
    fig.savefig('10to50_z_%.2fm' % pg.zplot[i] + '.png')    
    
    # plt.pause(1)
    if len(steps) > 1:
        if i<ppower.shape[0]-step:
            ax1.cla()
            ax12.cla()
            ax2.cla()
            ax3.cla()
            ax32.cla()
            ax4.cla()
            ax5.cla()

#%%%% animation for reports
import xkinterface.interface.PostGenesis13 as PostGenesis

workdir = f"/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4"  #"/Q_1nC_1Len_lscOFF/"
# os.chdir(workdir)
# workdir="/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_4nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_-20.0deg_phi1_0deg"
# workdir="/lustre/fs22/group/pitz/biaobin/2025/flattop_bcON_4nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_-16.0deg_phi1_0deg"
os.chdir(workdir)  
# # fname = "./02_gen4/gen4.seed5.nper5.fitdeg1.out.h5"
fname = "gen4.seed5.nper10.fitdeg3.out.h5"

# path = '/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4/all_lsc_off/'
path = '/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4/slsc_off_lscON/'
os.chdir(path)
fname = "gen4.out.h5"

# test for 170MeV beam
path = '/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/test_170MeV'

path='/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/Q_0.5nC_1Len_lscON'
os.chdir(path)
fname = "./gen4.seed5.nper10.fitdeg3.out.h5"

# workdir=f"/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4"
# os.chdir(workdir) 
# fname = "gen4.out.h5"

# tmp1="/afs/ifh.de/group/pitz/data/lixiangk/sim3/2023/SASEScan2/beam_200A_2.0nC_shotnoise_p/"
# os.chdir("/lustre/fs22/group/pitz/biaobin/xk_idealGaussFlattop")
# fname = tmp1+"pithz.1.out.h5"
# pg = PostGenesis(fname, debug=1, version=3, fig_ext = '.png')

pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

# fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(nrows = 5, figsize = (5, 14))
fig, [ax1, ax3, ax5] = plt.subplots(nrows = 3, figsize = (6, 10))
ax12 = ax1.twinx()
ax32 = ax3.twinx()

# plt.show()

name1 = 'intensity-farfield'
name2 = 'phase-farfield'

ppower = pg.get_fielddata('power')
p_mid = pg.get_fielddata(name1)
phi_mid = pg.get_fielddata(name2)
gamma = pg.get_beamdata('energy')
bunching = pg.get_beamdata('bunching')

xlim1 = 0
xlim2 = 100

step = 10
steps =np.arange(0, ppower.shape[0], step)
# steps = [40]
for i in steps:
    
    # ax1.plot(p_mid[i]/p_mid[i].max(), 'r-')
    # ax12.plot(ppower[i]/ppower[i].max(), 'b-')
    # ax1.plot(p_mid[i]/p_mid[i].max(), 'r-',label='intensity-farfield')
    # ax1.plot(ppower[i]/ppower[i].max(), 'b-',label='power')
    ax1.plot(ppower[i]/1e6, 'b-',label='power')
    
    # title = str.format('z = %.3f m,%0d undu.' % (pg.zplot[i], (pg.zplot[i]-0.105)/0.03))
    title = str.format('z = %.2f m' % (pg.zplot[i]))
    

    ax1.set_title(title)
    ax1.set_xlabel(r'# of slice')
    ax1.set_ylabel(r'On-axis power (MW)')
    #ax12.set_ylabel(r'Power (au)')
    # ax1.set_ylim(0, 2e7)
    ax1.grid()
    
    ax1.set_xlim(0, xlim2)
    # ax1.set_xlim(0,80)
    # ax1.set_xlim(5,60)
    # ax1.set_xlim(10,50)

    
    # ax1.set_ylim(0., 1.1)
    #ax1.set_yscale('log')
    ax1.legend()
    # ax12.set_ylim(0.0, 1.1)
    #ax12.set_yscale('log')
    
    # #=====================================
    # ax2.plot(gamma[0], 'r-',label='ini.')
    # ax2.plot(gamma[i], 'b-',label='step i')
    # ax2.set_xlabel('# of slice')
    # ax2.set_ylabel(r'$\gamma$')
    # # ax2.grid()
    # # ax2.set_ylim(32, 35)
    
    # ax2.set_xlim(0, xlim2)
    # # ax2.set_xlim(0, 125)
    # # ax2.set_xlim(10,50)
    # # ax2.set_xlim(0,30)
    # # ax2.set_xlim(20,100)
    # ax2.legend()
    
    #=====================================
    DE = (gamma[i]-gamma[0]) *pg.current*100e-6/2.998e8/1.6e-19 *0.511e6 *1.6e-19 *1e6 
    DE -= pg.cum_dE_lsc[i] *pg.current*100e-6/2.998e8/1.6e-19 #uJ/slice 
    
    # DE = (gamma[i]-gamma[0])*0.511e3 #keV
    #ax3.plot(DE/np.max(np.abs(DE)), 'r-')
    ax3.plot(DE, 'r-')
    
    ax3.set_xlim(0, xlim2)
    # ax3.set_xlim(10,50)
    # ax3.set_xlim(0,30)
    # ax3.set_xlim(20,100)
    
    # ax3.set_ylim(-2, 2)
    ax3.set_xlabel('# of slice')
    ax3.set_ylabel(r'Energy loss (uJ/slice)')
    ax3.grid(True)
    
    ax32.plot(pg.current, 'k-',label='current prof.')
    #ax32.set_ylabel(r'Current (A)')
    ax32.legend()
    
    # ax4.plot(bunching[i], 'r-')
    # ax4.set_xlabel('# of slice')
    # ax4.set_ylabel('Bunching')
    
    # ax4.set_ylim(1e-5, 1)
    # ax4.set_yscale('log')
    
    # ax4.set_xlim(0, xlim2)
    # # ax4.set_xlim(10,50)
    # # ax4.set_xlim(0,30)
    # # ax4.set_xlim(20,100)
    
    # ax4.grid(True)
    
    w, s = pg.get_spectrum(pg.zplot[i])
    #ax5.plot(w*1e6, s/np.max(s), '-')
    ax5.plot(w*1e6, s, '-')
    
    
    # w, s = calcSpectrum(p_mid[i], phi_mid[i])
    # ax5.plot(w*1e6, s/np.max(s), '-')
    
    # ax5.set_xlim(80, 120)
    #ax5.set_ylim(0, 1.1)
    # ax5.set_ylim(0, 2e13)
    
    ax5.set_xlabel(r'Wavelength ($\mu$m)')
    ax5.set_ylabel(r'Intensity (arb. unit)')
    ax5.grid(True)
    
    # display.display(plt.gcf())
    # display.clear_output(wait=True)

    fig.tight_layout()
    fig.savefig('report_z_%.2fm' % pg.zplot[i] + '.png')    
    
    # plt.pause(1)
    if len(steps) > 1:
        if i<ppower.shape[0]-step:
            ax1.cla()
            ax12.cla()
            # ax2.cla()
            ax3.cla()
            ax32.cla()
            # ax4.cla()
            ax5.cla()

#%%%% make gif video
path = '/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/'
path='/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4/slsc_off_lscON/'
path='/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4/all_lsc_off/'
os.chdir(path)

from PIL import Image
import glob

# os.remove("output.gif")
# 1. Load all PNG images in order
image_files = sorted(glob.glob("report_*.png"))  # Assumes images are named sequentially (e.g., frame1.png, frame2.png)

# image_files = image_files[::2]

# 2. Open each image and append to a list
images = []
for filename in image_files:
    img = Image.open(filename)
    images.append(img)

# 3. Save as GIF (adjust duration and loop as needed)
images[0].save(
    "output_3.gif",
    save_all=True,
    append_images=images[1:],  # Append remaining images
    duration=1000,  # Time per frame (milliseconds)
    loop=1,       # 0 = infinite loop, 1 = play once, etc.
)

print("GIF created successfully!")

#%%%% plots show seeding from tail
# 1nC, BC ON case
import xkinterface.interface.PostGenesis13 as PostGenesis

# for BC ON case
# ---------------
workdir = f"/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4"  #"/Q_1nC_1Len_lscOFF/"
os.chdir(workdir)  
fname = "gen4.seed5.nper10.fitdeg3.out.h5"

# for BC OFF case
#-----------------
workdir=f"/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/20241202-NonCompressor-ScanPhase/cores_016_impt_32320064_impz_32320064_Np_1.00e+07_phi2_-36.0deg/02_gen4"
os.chdir(workdir)
fname = "gen4.out.h5"

# workdir=f"/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4"
# os.chdir(workdir) 
# fname = "gen4.out.h5"

pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')


fig, [ax1, ax2] = plt.subplots(nrows= 2, figsize = (6, 8))
ax22 = ax2.twinx()

# plt.show()

name1 = 'intensity-farfield'
name2 = 'phase-farfield'

ppower = pg.get_fielddata('power')
p_mid = pg.get_fielddata(name1)
phi_mid = pg.get_fielddata(name2)
gamma = pg.get_beamdata('energy')
bunching = pg.get_beamdata('bunching')

xlim1 = 0
xlim2 = 60

step = 5
steps =np.arange(0, ppower.shape[0], step)
steps = [10]
for i in steps:

    ax1.plot(ppower[i]/1e6, 'b-')
    
    # title = str.format('z = %.2f m' % (pg.zplot[i]))
    
    period_i = (pg.zplot[i]-0.105)/0.03
    title = f"undu. period={period_i:.1f}"
    ax1.set_title(title)
    ax1.set_xlabel(r'# of slice')
    ax1.set_ylabel(r'THz power (MW)')
    #ax12.set_ylabel(r'Power (au)')
    # ax1.set_ylim(0, 2e7)
    ax1.grid()
    
    ax1.set_xlim(0, xlim2)
    ax1.legend()    

    DE = (gamma[i]-gamma[0]) *pg.current*100e-6/2.998e8/1.6e-19 *0.511e6 *1.6e-19 *1e6 
    DE -= pg.cum_dE_lsc[i] *pg.current*100e-6/2.998e8/1.6e-19 #uJ/slice 
    # DE = (gamma[i]-gamma[0])*0.511e3 #keV
    #ax3.plot(DE/np.max(np.abs(DE)), 'r-')
    ax2.plot(DE, 'r-') #,label="Energy loss")
    
    ax2.set_xlim(0, xlim2)
    ax2.set_xlabel('# of slice')
    ax2.set_ylabel(r'Energy loss from THz (uJ/slice)')
    ax2.grid(True)
    ax2.legend()
    
    ax22.plot(pg.current, 'k-',label="current profile")
    ax22.set_ylabel(r'current (A)')
    ax22.legend()
    
    fig.tight_layout()
    fig.savefig('10to50_z_%.2fm_report' % pg.zplot[i] + '.png')    
    
#%%%% one plot
# 1nC, BC ON case
import xkinterface.interface.PostGenesis13 as PostGenesis

# for BC ON case
# ---------------
workdir = f"/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4"  #"/Q_1nC_1Len_lscOFF/"
os.chdir(workdir)  
fname = "gen4.seed5.nper10.fitdeg3.out.h5"

# for BC OFF case
#-----------------
# workdir=f"/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/20241202-NonCompressor-ScanPhase/cores_016_impt_32320064_impz_32320064_Np_1.00e+07_phi2_-36.0deg/02_gen4"
# os.chdir(workdir)
# fname = "gen4.out.h5"

# workdir=f"/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4"
# os.chdir(workdir) 
# fname = "gen4.out.h5"

pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')


fig, ax1 = plt.subplots(nrows= 1, figsize = (6, 5))
ax12 = ax1.twinx()

# plt.show()

name1 = 'intensity-farfield'
name2 = 'phase-farfield'

ppower = pg.get_fielddata('power')
p_mid = pg.get_fielddata(name1)
phi_mid = pg.get_fielddata(name2)
gamma = pg.get_beamdata('energy')
bunching = pg.get_beamdata('bunching')

xlim1 = 0
xlim2 = 60

step = 5
steps =np.arange(0, ppower.shape[0], step)
steps = [10]
for i in steps:
    period_i = (pg.zplot[i]-0.105)/0.03
    title = f"undu. period={period_i:.1f}"
    ax1.set_title(title)

    DE = (gamma[i]-gamma[0]) *pg.current*100e-6/2.998e8/1.6e-19 *0.511e6 *1.6e-19 *1e6 
    DE -= pg.cum_dE_lsc[i] *pg.current*100e-6/2.998e8/1.6e-19 #uJ/slice 
    # DE = (gamma[i]-gamma[0])*0.511e3 #keV
    #ax3.plot(DE/np.max(np.abs(DE)), 'r-')
    ax1.plot(DE*1e3, 'r-') #,label="Energy loss")
    
    ax1.set_xlim(0, xlim2)
    ax1.set_xlabel('# of slice')
    ax1.set_ylabel(r'Energy loss from THz (nJ/slice)')
    ax1.grid(True)
    ax1.legend()
    
    ax12.plot(pg.current, 'k-',label="current profile")
    ax12.set_ylabel(r'current (A)')
    ax12.legend()
    
    fig.tight_layout()
    fig.savefig('10to50_z_%.2fm_report' % pg.zplot[i] + '.png')     

#%%% peak power - phi2
# for flattop beam
bc_statl = ["OFF","ON"]
phi2l = np.arange(-6,-32-1,-2)

phi2l = np.arange(0,-26-1,-2)

plt.figure(figsize=(8,6))
for bc_stat in bc_statl:
    peakpower = []
    for phi2 in phi2l:
        workdir= f"/lustre/fs22/group/pitz/biaobin/2025/flattop_bc{bc_stat}_4nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_{phi2}.0deg_phi1_0deg"
        workdir =f"/lustre/fs22/group/pitz/biaobin/2025/flattop_bc{bc_stat}_4nC/Np_5.24e+06_phi1_0deg_phi2_{phi2}deg"
        workdir =f"/lustre/fs22/group/pitz/biaobin/2025/flattop_bc{bc_stat}_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{phi2}.0deg"

        os.chdir(workdir)   
        
        if bc_stat=="OFF":
            fname = "./02_gen4/gen4.seed5.nper5.fitdeg1.out.h5"
        else:
            fname = "./02_gen4/gen4.seed5.nper10.fitdeg3.out.h5"
            
        pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')
        
        tmp = np.max( pg.get_tpower(at=pg.Lsat) )
        peakpower.append(  tmp )
        
    peakpower = np.array(peakpower)
    plt.plot(phi2l, peakpower/1e6, label=f"bc {bc_stat}")

plt.legend()
plt.grid()
plt.xlabel("phi2 (deg)")
plt.ylabel("THz peak power (MW)")

#%%%% power@exit - phi2
# for flattop beam
bc_statl = ["OFF","ON"]
phi2l = np.arange(-6,-32-1,-2)

phi2l = np.arange(0,-26-1,-2)

plt.figure(figsize=(8,6))
for bc_stat in bc_statl:
    peakpower = []
    for phi2 in phi2l:
        workdir= f"/lustre/fs22/group/pitz/biaobin/2025/flattop_bc{bc_stat}_4nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_{phi2}.0deg_phi1_0deg"
        workdir =f"/lustre/fs22/group/pitz/biaobin/2025/flattop_bc{bc_stat}_4nC/Np_5.24e+06_phi1_0deg_phi2_{phi2}deg"
        workdir =f"/lustre/fs22/group/pitz/biaobin/2025/flattop_bc{bc_stat}_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{phi2}.0deg"

        os.chdir(workdir)   
        
        if bc_stat=="OFF":
            fname = "./02_gen4/gen4.seed5.nper5.fitdeg1.out.h5"
        else:
            fname = "./02_gen4/gen4.seed5.nper10.fitdeg3.out.h5"
            
        pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')
        
        tmp = pg.zenergy[-1]
        peakpower.append(  tmp )
        
    peakpower = np.array(peakpower)
    plt.plot(phi2l, peakpower*1e3, label=f"bc {bc_stat}")

plt.legend()
plt.grid()
plt.xlabel("phi2 (deg)")
plt.ylabel("THz energy (mJ)")


#%%%% peak power of gaussian beam
bc_statl = ["ON", "OFF"]
phi2l = np.arange(-16,-40-1,-2)

plt.figure(figsize=(8,6))
for bc_stat in bc_statl:
    peakpower = []
    for phi2 in phi2l:
        #
        workdir= f"/lustre/fs22/group/pitz/biaobin/2025/gauss_bc{bc_stat}_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{phi2}.0deg"
        os.chdir(workdir)   
        
        if bc_stat=="OFF":
            fname = "./02_gen4/gen4.seed5.nper5.fitdeg1.out.h5"
        else:
            fname = "./02_gen4/gen4.seed5.nper10.fitdeg3.out.h5"
            
        pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')
        
        tmp = np.max( pg.get_tpower(at=pg.Lsat) )
        peakpower.append(  tmp )
        
    peakpower = np.array(peakpower)
    plt.plot(phi2l, peakpower/1e6, label=f"bc {bc_stat}")

plt.legend()
plt.grid()
plt.xlabel("phi2 (deg)")
plt.ylabel("THz peak power (MW)")

#%%% LSC along undulator
workdir= f"/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg"

os.chdir(workdir) 

fname = "./02_gen4/gen4.seed5.nper10.fitdeg3.out.h5"
# fname = "./02_gen4/gen4.seed5.nper5.fitdeg1.out.h5"

# pg1 = PostGenesis(fname, debug=1, version=4, fig_ext = '.png')
# pg1.plot_tpower(at=2.6)
# pg1.plot_spectrum(at=2.6)
# pg1.plot_current()

ff = h5py.File(fname,'r')

# plt.figure(figsize=(8,6))
# plt.plot(ff["Beam"]["LSCfield"][240])


lsc = ff["Beam"]["LSCfield"]
current = ff["Beam"]["current"][:].flatten()

nstep = lsc.shape[0]

dgam_lsc =[]
for j in np.arange(nstep):

    lsc_s = lsc[j].flatten()
    dgam_lsc.append( lsc_s*1.6e-19 *current*100e-6/2.998e8 *0.015 *1e6) #[uJ/slice]
    
    dgam_lsc.append(  )

dgam_lsc = np.array(dgam_lsc)

cum_dgam_lsc = np.cumsum(dgam_lsc, axis=0)


plt.figure()

plt.plot(cum_dgam_lsc[240],'-')

#%%% beam phase evo
from xkinterface.interface import Plot as xkplt
import impzpy.post as pzplot

# xkplt.plot2d(pha.zt,pha.dgam,xlabel="z (ps)",ylabel=r"$\Delta \gamma$")

# # workdir=f"/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4/slsc_off_lscON"
# workdir=f"/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4/slsc_off"
# # workdir=f"/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4"
# # workdir=f"/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/20241202-NonCompressor-ScanPhase/cores_016_impt_32320064_impz_32320064_Np_1.00e+07_phi2_-32.0deg/02_gen4"
# os.chdir(workdir)

# path = '/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/test_170MeV'
# # os.chdir(path)

# path='/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/Q_0.5nC_1Len_lscON'
# path='/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/Q_1nC_2Len_lscON/'
# path="/lustre/fs22/group/pitz/biaobin/2025/flattop_bcOFF_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-6.0deg/02_gen4/addPeak"

# path="/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/test_1THz"
# os.chdir(path)

import h5py

sampleFreq=20

step = [int(0.1/0.15)]
# step = [int(3.0/0.15)]
# step = [int(2.25/0.15)]

# for jj in range(25):
for jj in [10]:
# for jj in step:
    fname = 'gen4.seed5.nper10.fitdeg3.'+str(jj*10)+'.par.h5'
    # fname = 'gen4.'+str(jj*10)+'.par.h5'

    f = h5py.File(fname,'r')

    position_s = jj*10*0.015 #m
    s1 = 0
    s2 = 45
    # s2 = 60 #128
    # s2 = 90
    # s2 = 100
    
    dlamda = s1
    lamdas = 1 
    x= np.array([])
    y= np.array([])
    # for j in range(90):
    for j in np.arange(s1,s2,1):
        if j+1 < 10:
            name = 'slice000'+'00'+str(j+1)
        elif j+1 < 100:
            name = 'slice000'+'0'+str(j+1)
        elif j+1 < 1000:
            name = 'slice000'+str(j+1)
    
        x = np.append(x,f[name]['theta'][:]/(2*np.pi)*lamdas+dlamda)
        y = np.append(y,f[name]['gamma'][:])
        # plt.plot(f[name]['theta'][:]/np.pi*180+dphi,f[name]['gamma'][:],'.r')
        dlamda = dlamda +lamdas
    
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
    pzplot.plot_dotpha(x[::sampleFreq],y[::sampleFreq],bins=256*2*2)
    plt.xlim([20,30])
    plt.xlim([22,40])
    plt.xlim([40,75])
    plt.xlim([5,44])
    plt.xlabel(r'$z/\lambda_s$')
    plt.ylabel(r'$\gamma$')
    plt.grid()
    

    plt.savefig(f"{position_s:.2f}m_beamPhase.png")

    # plt.close()

    # plt.subplot(1,2,1)
    # xkplt.plot2d(x, y,bins=[1000,1000])
#%%%% make gif video
# path = '/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4/slsc_off_lscON/'
path = '/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/02_gen4/all_lsc_off/'
os.chdir(path)

from PIL import Image
import glob

# os.remove("output.gif")
# 1. Load all PNG images in order
image_files = sorted(glob.glob("*_beamPhase.png"))  # Assumes images are named sequentially (e.g., frame1.png, frame2.png)
gif_name = "phase_evo.gif"

# 2. Open each image and append to a list
images = []
for filename in image_files:
    img = Image.open(filename)
    images.append(img)

# 3. Save as GIF (adjust duration and loop as needed)
images[0].save(
    gif_name,
    save_all=True,
    append_images=images[1:],  # Append remaining images
    duration=1000,  # Time per frame (milliseconds)
    loop=1,       # 0 = infinite loop, 1 = play once, etc.
)

print("GIF created successfully!")












    
