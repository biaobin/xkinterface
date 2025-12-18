import os
import numpy as np
import matplotlib.pyplot as plt
from impzpy import post
from scipy.interpolate import interp1d
from cycler import cycler
import h5py
import xkinterface.interface.PostGenesis13 as PostGenesis
import sys
from startup import *


# Define markers, colors, and line styles
markers = ['o', 's', 'p', 'v', 'd', '*', '+', 'D', '>', 'h', 'H', '^', '<', '.', ',']
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '#1f77b4', '#ff7f0e', '#2ca02c', 
          '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22']
linestyles = ['-', '--', '-.', ':', (0, (1, 1)), (0, (5, 1)), (0, (3, 1, 1, 1)), 
              (0, (3, 5, 1, 5)), (0, (2, 2, 3, 3, 5, 5)), (0, (1, 10)), (0, (5, 10))]
# Ensure all lists are the same length by repeating or trimming
num_styles = min(len(markers), len(colors), len(linestyles))
markers = markers[:num_styles]
colors = colors[:num_styles]
linestyles = (linestyles * (num_styles // len(linestyles) + 1))[:num_styles]  # Repeat & trim
# Set the default property cycle
plt.rc('axes', prop_cycle=cycler(marker=markers) + cycler(color=colors) + cycler(linestyle=linestyles))


# basedir = '/afs/ifh.de/group/pitz/data/biaobin/sim1/2025/00_sase_chicane/Q_1_1.5_2nC/bcON_1nC'
basedir = '/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC'
basedir = '/mnt/f/simu_2025/test_17MeV_one4one_scOFF'
os.chdir(basedir)

#%% impt section
#%%% longi-phase & current profile
basedir = r'/lustre/fs22/group/pitz/biaobin/2025/twobeam_bcON_1nC_delay3.0sigt_Q2_100pC_sigt2_250fs'
# basedir = r'/lustre/fs22/group/pitz/biaobin/2025/gauss_bcOFF_1nC'
os.chdir(basedir)

path=r'cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_-32.0deg/00_impt'
os.chdir(path)

id = 90
plt.figure(figsize=(10,4))

fname1 = f"fort.{id}"
fname2 = f"fort.{id+100}"

fort1001= np.loadtxt(fname1)
fort11001 =np.loadtxt(fname2)

plt.subplot(1,2,1)
# plt.hist2d(fort1001[:,4]*1e3,fort1001[:,5],bins=200,cmap='jet')
plt.plot(fort1001[:,4]*1e3,fort1001[:,5],'.',markersize=1)
plt.xlabel('z (mm)')
plt.ylabel("gambetz")

plt.subplot(1,2,2)
plt.plot(fort11001[:,0]*1e3,fort11001[:,2],'-r')
plt.xlabel('z (mm)')
plt.ylabel('current (A)')


#%%% movie
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

# Define the directory path
bc_stat= "OFF"
path = r'/afs/ifh.de/group/pitz/data/biaobin/sim2/2025/twobeam_bcOFF_1nC/cores_016_impt_32320064_impz_32320064_Np_1.31e+06_phi2_0.0deg/00_impt'
path=r'/mnt/d/NextCloud/subjects/2025/01_twoBeam/twobeam_bcOFF_1nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_0.0deg/00_impt'

path = f'/lustre/fs22/group/pitz/biaobin/2025/flattop_bc{bc_stat}_1nC'

path = f'/lustre/fs22/group/pitz/biaobin/2025/flattop_bc{bc_stat}_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-16.0deg/00_impt'

path='/lustre/fs22/group/pitz/biaobin/test_2025_0216/twobeam_bcON_1nC_delay2.0sigt_Q2=100pC_sigt2=500fs/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_-20.0deg/00_impt'
path='/lustre/fs22/group/pitz/biaobin/test_2025_0216/twobeam_bcON_1nC_delay2.0sigt_Q2100pC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_-20.0deg/00_impt'
path='/lustre/fs22/group/pitz/biaobin/test_2025_0216/twobeam_bcOFF_1nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_-20.0deg/00_impt'
os.chdir(path)

# Define idd range
idd = np.arange(1000, 1496, 1)

# Create figure and subplots
fig, ax = plt.subplots(1, 2, figsize=(10, 4))

# Subplot 1: Scatter plot
scatter, = ax[0].plot([], [], '.', label="Phase Space")
ax[0].set_xlabel('z (mm)')
ax[0].set_ylabel("gambetz")

# Subplot 2: Line plot
line, = ax[1].plot([], [], '-r', label="Current Profile")
ax[1].set_xlabel('z (mm)')
ax[1].set_ylabel('current (A)')

# Update function for animation
def update(frame):
    id = idd[frame]  # Current id value
    fname1 = f"fort.{id}"
    fname2 = f"fort.{id+10000}"

    # Check if files exist
    if not os.path.exists(fname1) or not os.path.exists(fname2):
        print(f"Skipping {id}, files not found.")
        return scatter, line

    # Load data
    fort1001 = np.loadtxt(fname1)
    fort11001 = np.loadtxt(fname2)

    # Update scatter plot
    scatter.set_data(fort1001[:, 4] * 1e3, fort1001[:, 5])
    ax[0].set_xlim(np.min(fort1001[:, 4]) * 1e3, np.max(fort1001[:, 4]) * 1e3)
    ax[0].set_ylim(np.min(fort1001[:, 5]), np.max(fort1001[:, 5]))

    # Update line plot
    line.set_data(fort11001[:, 0] * 1e3, fort11001[:, 2])
    ax[1].set_xlim(np.min(fort11001[:, 0]) * 1e3, np.max(fort11001[:, 0]) * 1e3)
    ax[1].set_ylim(np.min(fort11001[:, 2]), np.max(fort11001[:, 2]))

    return scatter, line

# Create animation
ani = animation.FuncAnimation(fig, update, frames=len(idd), interval=200, blit=False)
plt.show()
# ani.save("animation.mp4", fps=10, dpi=150)

#%% impz section

#%%% phi2-current
# profile="gauss"
profile="flattop"

bc_stat="ON"
charge=2

if profile =="gauss" or profile=="flattop":
    core=6464
    Np="5.24e+06"
# elif profile =="flattop":
#     core=3232
#     Np="6.55e+05"

basedir = f'/lustre/fs22/group/pitz/biaobin/2025/{profile}_bc{bc_stat}_{charge}nC'
os.chdir(basedir)

tmp1 = f'cores_016_impt_{core}0064_impz_{core}0064_Np_{Np}_phi2_' 
tmp2 = 'deg/01_impz/fort.11003'
           
degl = np.arange(0, -26-1, -2)
# degl = np.arange(-16, -40-1, -2)


plt.figure(figsize=(16,6))

# dat = np.loadtxt('/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/202411-ScanBooster-S2E/20241130_grid_3232064/phi2_-32deg/01_impz/fort.11000')
# plt.plot(dat[:,0]*1e3,dat[:,2],'-',label='cds exit')
plt.subplot(1,2,1)
peak_cur = []

j=0
for jj in degl:
    tmp = tmp1 +"{:.1f}".format(jj) +tmp2
    dat = np.loadtxt(tmp)
    peak_cur.append(np.max(dat[:,2]))
    plt.plot(dat[:,0]/2.998e8*1e12,dat[:,2],'-',label=str(jj)+'deg')
    j +=1

plt.title(f"BC {bc_stat}, {charge} nC")    
plt.xlabel('z (ps)')
plt.ylabel('current @und. ent. (A)')
plt.legend()
plt.grid()
# plt.xlim([-10,10])

plt.subplot(1,2,2)
plt.plot(degl,peak_cur)
plt.xlabel('phi2 (deg)')
plt.ylabel('peak current @und. ent. (A)')
plt.grid()


dat = np.column_stack((degl,peak_cur))
np.savetxt("degl_peakcurrent.txt", dat)
    

plt.savefig("phi2-current.png")

#phi2-core_current
#-----------------
# get the core-current manually

# for jj in degl:
#     plt.figure(figsize=(8,6))
#     tmp = tmp1 +"{:.1f}".format(jj) +tmp2
#     dat = np.loadtxt(tmp)
    
#     plt.plot(dat[:,0]/2.998e8*1e12,dat[:,2],'-',label=str(jj)+'deg')
#     plt.title("phi2="+str(jj))
    
# plt.xlabel('z (ps)')
# plt.ylabel('current @und. ent. (A)')
# plt.legend()

# get the core-current manually
# core_current=[285, 660, 1174, 863, 325, 238, 166, 129, 108, 97, 83, 77, 71]
# core_current.reverse()

# plt.figure(figsize=(8,6))
# plt.plot(degl,core_current,'o-')
# plt.xlabel('phi2 (deg)')
# plt.ylabel('core sec. current @und. ent. (A)')
# plt.show()    
# plt.savefig("phi2--core-current.png")

#%%%% plot together, different charge, phi2-peakCurrent
rootdir = "/lustre/fs22/group/pitz/biaobin/2025/"
profile="flattop"

# bc_statl=["OFF","ON"]
bc_statl=["ON"]

chargel = [1,1.5,2]

if profile =="gauss" or profile=="flattop":
    core=6464
    Np="5.24e+06"
# elif profile =="flattop":
#     core=3232
#     Np="6.55e+05"

plt.figure(figsize=(8,6))

for bc_stat in bc_statl:
    for charge in chargel:
        basedir = f'/lustre/fs22/group/pitz/biaobin/2025/{profile}_bc{bc_stat}_{charge}nC/'
        
        fname = basedir + "degl_peakcurrent.txt"
        dat = np.loadtxt(fname)
        
        
        plt.plot(dat[:,0],dat[:,1],label=f"BC {bc_stat:>3}, {charge:.2f} nC")
        
plt.legend()        
plt.title(f"{profile} profile")        
plt.xlabel("phi2 (deg)")
plt.ylabel("peak current (A)")
plt.grid("ON")

plt.savefig(f"{rootdir}{profile}_{bc_statl}_phi2_peakcurrent.png")


#%%% phi2-sigz,core-current,sigE,enxeny
# profile="flattop"

profile="gauss"
bc_stat = "ON"
charge= 2

basedir = f'/lustre/fs22/group/pitz/biaobin/2025/{profile}_bc{bc_stat}_{charge}nC'
# basedir = f'/afs/ifh.de/group/pitz/data/biaobin/sim2/2025_bugPhaseV/bcOFF_1nC/'
os.chdir(basedir)

if profile =="gauss":  # or profile=="flattop":
    core=6464
    Np="5.24e+06"
    degl = np.arange(-16,-41,-2)

elif profile =="flattop":
    core=3232
    Np="6.55e+05"
    degl = np.arange(0, -40 -1, -2)
    
# degl = [-36]

sigzl = []
enxl = []
enyl = []
sigEl = []
peakCur=[]
for deg in degl:    
    tmppath =f"cores_016_impt_{core}0064_impz_{core}0064_Np_{Np}_phi2_{deg:.1f}deg/01_impz/"
    
    tmp = post.twissplot(path=tmppath)
    # plt.figure(figsize=(6,5))
    # plt.plot(tmp.s,tmp.sigz,'r-')
    
    interp_func_sigz= interp1d(tmp.s, tmp.sigz, kind='linear', fill_value="extrapolate")
    interp_func_enx = interp1d(tmp.s, tmp.enx, kind='linear', fill_value="extrapolate")
    interp_func_eny = interp1d(tmp.s, tmp.eny, kind='linear', fill_value="extrapolate")
    interp_func_sigE = interp1d(tmp.s, tmp.sigE, kind='linear', fill_value="extrapolate")
    
    s_undu = 22.94
    sigzl.append(interp_func_sigz(s_undu))
    enxl.append(interp_func_enx(s_undu))
    enyl.append(interp_func_eny(s_undu))
    sigEl.append(interp_func_sigE(s_undu))
    
sigzl = np.array(sigzl)*1e-3/2.998e8*1e12  #ps        
plt.figure(figsize=(24,6))
plt.subplot(1,3,1)
plt.plot(degl,sigzl,'o-')    
plt.xlabel('phi2 (deg)')
plt.ylabel("sigz (ps)")

# plt.subplot(1,4,2)
# get the core-current manually
# core_current=[285, 660, 1174, 863, 325, 238, 166, 129, 108, 97, 83, 77, 71]
# core_current.reverse()
# plt.plot(degl,core_current,'o-')
# plt.xlabel('phi2 (deg)')
# plt.ylabel('core sec. current @und. ent. (A)')

plt.subplot(1,3,2)
plt.plot(degl,np.array(sigEl)*1e3,'o-')    
plt.xlabel('phi2 (deg)')
plt.ylabel("sigE (KeV)")

plt.subplot(1,3,3)
plt.plot(degl,enxl,'o-',label='enx')    
plt.xlabel('phi2 (deg)')
plt.ylabel("enx (mm mrad)")

plt.plot(degl,enyl,'^--',label='eny')    
plt.xlabel('phi2 (deg)')
plt.ylabel("norm. emit. (mm mrad)")
plt.legend()
plt.savefig("phi2-sigz-sigE-enxeny.png")


#%%% phi2-longiPhase
bc_stat = "OFF"
profile="flattop"

charge= 1

basedir = f'/lustre/fs22/group/pitz/biaobin/2025/test_flattop_rs_3ps_22psTop/{profile}_bc{bc_stat}_{charge}nC'
# basedir = f'/lustre/fs22/group/pitz/biaobin/2025/{profile}_bc{bc_stat}_{charge}nC_delay3.0sigt_Q2_100pC_sigt2_250fs'
# basedir = f'/afs/ifh.de/group/pitz/data/biaobin/sim2/2025_bugPhaseV/bcOFF_1nC/'

os.chdir(basedir)

if profile =="gauss" or profile=="flattop":
    core=6464
    Np="5.24e+06"
elif profile=="twobeam":
    core=3232
    Np="6.55e+05"

# degl = np.arange(-16, -40 -1, -2)
degl = [-28,-30,-32,-34]

degl = [-14]
for deg in degl:    
    #---------------------------------------------------------------
    
    tmp = f"cores_016_impt_{core}0064_impz_{core}0064_Np_{Np}_phi2_{deg}.0deg/01_impz/"
    # tmp = f"cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{deg}.0deg/01_impz/"
    
    # rms-evo
    twi = post.twissplot(path=tmp)
    
    twi.plot_rms(enxy='ON')
    plt.savefig("phi2="+str(deg)+"deg-"+"rmsevo.png")
    
    # continue

    #jet plot
    #---------------------------------------------------------------
    def zdgamplot(fname1):
        # fname1 = 'fort.1000'
        # fname2 = 'fort.11000'
        fname2 = fname1.split('.')[0]+'.'+str(10000+int(fname1.split('.')[1]))
        
        fort1001= np.loadtxt(tmp+fname1,skiprows=1)
        fort11001 =np.loadtxt(tmp+fname2)
        
        
        plt.figure(figsize=(12,5))
        
        plt.subplot(1,2,1)
        plt.title(f"BC {bc_stat}, phi2="+str(deg)+" deg, "+fname1)
        plt.hist2d(fort1001[:,4]/2.998e8*1e12,fort1001[:,5],bins=300,cmap='jet')
        plt.xlabel('z (ps)')
        plt.ylabel(r'$\Delta\gamma$')
        
        plt.subplot(1,2,2)
        plt.plot(fort11001[:,0]/2.998e8*1e12,fort11001[:,2],'-r')
        plt.xlabel('z (ps)')
        plt.ylabel('current (A)')
        
        plt.savefig("phi2="+str(deg)+"deg-"+fname1+"-zdgam.png")
    
    fname1 = 'fort.1000'
    zdgamplot(fname1)
    
    if bc_stat == "ON":
        fname1 = 'fort.1001'
        zdgamplot(fname1)
        
        fname1 = 'fort.1002'
        zdgamplot(fname1)
    
    fname1 = 'fort.1003'
    zdgamplot(fname1)
    
    # dot plot
    #---------------------------------------------------------------
    #plot together
    def zdgam_dotplot(fname1):
        # fname1 = 'fort.1000'
        # fname2 = 'fort.11000'
        fname2 = fname1.split('.')[0]+'.'+str(10000+int(fname1.split('.')[1]))
        
        fort1001= np.loadtxt(tmp+fname1,skiprows=1)
        fort11001 =np.loadtxt(tmp+fname2)
        
        plt.subplot(1,4,1)
        plt.title(f"BC {bc_stat}, phi2="+str(deg)+"deg")
        plt.plot(fort1001[::10,4]/2.998e8*1e12,fort1001[::10,5],'.')
        plt.xlabel('z (ps)')
        plt.ylabel(r'$\Delta\gamma$')
        
        plt.subplot(1,4,2)
        plt.plot(fort11001[:,0]/2.998e8*1e12,fort11001[:,2],'-')
        plt.xlabel('z (ps)')
        plt.ylabel('current (A)')
        
        if fname1 in ["fort.1001","fort.1002"]:
            if fname1 == "fort.1001":
                labl = "before BC"
            elif fname1 == "fort.1002":
                labl = "after BC"

            plt.subplot(1,4,3)
            plt.plot(fort11001[:,0]/2.998e8*1e12,fort11001[:,3]*1e6,'-',label=labl+': enx')
            plt.plot(fort11001[:,0]/2.998e8*1e12,fort11001[:,4]*1e6,'--',label=labl+': eny')
            plt.xlabel('z (ps)')
            plt.ylabel('slice enx & eny (um rad)') 
            plt.legend()
            
            plt.subplot(1,4,4)
            plt.plot(fort11001[:,0]/2.998e8*1e12,fort11001[:,6]/1e3,'-',label=labl)
            plt.xlabel('z (ps)')
            plt.ylabel('slice energy spread (keV)')
            plt.legend()    # os.chdir(tmp)

    
    plt.figure(figsize=(28,5))
    
    zdgam_dotplot('fort.1000')
    if bc_stat == "ON":
        zdgam_dotplot('fort.1001')
        zdgam_dotplot('fort.1002')
    zdgam_dotplot('fort.1003')
    
    plt.savefig("phi2="+str(deg)+"deg-"+"zdgam_0123.png",bbox_inches='tight')
    
    plt.close("all")


    
#%%%% single plot 
deg = -32

bc_stat = "ON"
profile="gauss"
charge= 1

basedir = f'/lustre/fs22/group/pitz/biaobin/2025/{profile}_bc{bc_stat}_{charge}nC'
os.chdir(basedir)
  
tmp1 = 'cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_' 
tmp2 = 'deg/01_impz/'
tmp3 = tmp1 +"{:.1f}".format(deg)+"deg"
tmp =tmp3+"/01_impz/"

def plot_zdgam_jet_dot(fname1):
    fort1001 = np.loadtxt(tmp+fname1,skiprows=1)
    plt.figure(figsize=(15,6))
    
    plt.subplot(1,2,1)
    plt.title("phi2="+str(deg)+" deg, "+fname1)
    plt.hist2d(fort1001[:,4]*1e3,fort1001[:,5],bins=300,cmap='jet')
    plt.xlabel('z (mm)')
    plt.ylabel(r'$\Delta\gamma$')
    
    plt.subplot(1,2,2)
    plt.plot(fort1001[:,4]*1e3,fort1001[:,5],'.')
    plt.xlabel('z (mm)')
    plt.ylabel(r'$\Delta\gamma$')    
    
    plt.savefig("phi2="+str(deg)+"deg-"+fname1+"-zdgam-jet-dot.png")

fname1 = "fort.1001"    
plot_zdgam_jet_dot(fname1)    
    
fname1 = "fort.1002"    
plot_zdgam_jet_dot(fname1) 

fname1 = "fort.1003"    
plot_zdgam_jet_dot(fname1) 

#%% genesis4 section
#%%% (method-1) phi2-gainCurve-Energy
#%%%% given gen4 path
workdir ='/lustre/fs22/group/pitz/biaobin/2025/test_flattop_rs_3ps_22psTop/flattop_bcOFF_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-8.0deg/02_gen4/copenChirp'
os.chdir(workdir)

fname = "gen4.seed5.nper5.fitdeg1.out.h5"
pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

fname2 = "../gen4.seed5.nper5.fitdeg1.out.h5"
pg2 = PostGenesis(fname2, debug=0, version=4, fig_ext = '.png')

plt.figure(figsize=(6,5))
plt.plot(pg.zplot, pg.zenergy*1e6,'-',label='compensated')
plt.yscale("log")

plt.plot(pg2.zplot, pg2.zenergy*1e6,'--',label='No compensation')
plt.yscale("log")
plt.grid()
plt.xlabel("s (m)")
plt.ylabel("THz energy (uJ)")
plt.legend()


#%%%% single case 
#read a single h5 output to plot the result
rootdir = "/lustre/fs22/group/pitz/biaobin/2025/"

# profile = "gauss"
profile = "flattop"

charge = 2
bc_stat = "OFF"

basedir=f'/lustre/fs22/group/pitz/biaobin/2025/{profile}_bc{bc_stat}_{charge}nC'
os.chdir(basedir)

# degl = np.arange(-16, -40 -1, -2)
degl = np.arange(0, -26 -1, -2)
# degl = [-8]

plt.figure(figsize=(8,6))

maxE = []
j=0
for deg in degl:  
    tmp = "cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{:.1f}deg/02_gen4/".format(deg)
    # tmp = "cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_{:.1f}deg/02_gen4/".format(deg)
    
    if bc_stat == "ON":
        fname = "gen4.seed5.nper10.fitdeg3.out.h5"
    else:
        fname = "gen4.seed5.nper5.fitdeg1.out.h5"
    pg = PostGenesis(tmp+fname, debug=0, version=4, fig_ext = '.png')

    maxE.append(np.max(pg.zenergy))

    # plt.plot(pg.zplot[::5],pg.zenergy[::5]*1e6,'-',color=colors[j],marker=markers[j],label="phi2="+str(deg)+"deg")
    plt.plot(pg.zplot[::5],pg.zenergy[::5]*1e6,label="phi2="+str(deg)+"deg",markersize=8)
    
    plt.yscale('log')
    
    plt.legend()
    plt.xlabel('z (m)')
    plt.ylabel('Energy (uJ)')    
    plt.title("")
    j +=1
plt.grid()
plt.title(f"{profile}, BC {bc_stat}, {charge:.2f} nC")
# plt.savefig(f"{rootdir}phi2-allGainCurves-Energy"+fname+".png")
plt.savefig(f"phi2-allGainCurves-Energy"+fname+".png")



# plt.figure(figsize=(8,6))
# plt.plot(degl,np.array(maxE)*1e6,'o-')
# plt.xlabel("phi2 (deg)")
# plt.ylabel("THz energy@undu exit (uJ)")
# plt.grid()
# plt.savefig("phi2-THzEnergy-"+fname+".png")

#%%% phi2 - power/energy for several charges
rootdir = "/lustre/fs22/group/pitz/biaobin/2025/"

# profilel = ["gauss","flattop"]
profilel = ["flattop"]
# bc_statl = ["OFF","ON"]
bc_statl = ["ON"]

chargel = [1,1.5,2]
degl = np.arange(0, -26 -1, -2)
seed= 5

keyword="power"

maxE = []
for profile in profilel:
    for bc_stat in bc_statl:
        for charge in chargel:
            
            basedir=f'/lustre/fs22/group/pitz/biaobin/2025/{profile}_bc{bc_stat}_{charge}nC'
            os.chdir(basedir)
            
            if profile =="gauss" or profile=="flattop":
                core=6464
                Np="5.24e+06"
            # elif profile =="flattop":
            #     core=3232
            #     Np="6.55e+05"
                
            tmp_maxE = []    
            for deg in degl:        
                tmp = f"cores_016_impt_{core}0064_impz_{core}0064_Np_{Np}_phi2_{deg:.1f}deg/02_gen4/"
                
                if bc_stat == "OFF":
                    fname = f"gen4.seed{seed}.nper5.fitdeg1.out.h5"
                elif bc_stat == "ON":
                    fname = f"gen4.seed{seed}.nper10.fitdeg3.out.h5"
                # print("fname=",fname)
                pg = PostGenesis(tmp+fname, debug=0, version=4, fig_ext = '.png')
                if keyword=="energy":
                    tmp_maxE.append(np.max(pg.zenergy)*1e6 ) #uJ
                elif keyword == "power":
                    tmp_maxE.append(np.max(pg.zpower)/1e6)   #MW
            maxE.append(tmp_maxE)
#======            
plt.figure(figsize=(8*1.2,6*1.2))
j=0   
for profile in profilel:
    for bc_stat in bc_statl:
        for charge in chargel:     
            maxE = np.array(maxE)                
            plt.plot(degl,maxE[j],label=f"{profile}, BC {bc_stat:>3}, {charge:.1f} nC")
            plt.xlabel("phi2 (deg)")
            if keyword == "energy":
                plt.ylabel(f"THz {keyword} (uJ)")
            elif keyword == "power":
                plt.ylabel(f"THz {keyword} (MW)")
            plt.legend()
            j += 1
plt.grid()

os.chdir(rootdir)
plt.savefig(f"{rootdir}{'_'.join(profilel)}_{'_'.join(bc_statl)}_phi2_THz{keyword}_seed{seed}.png")

#%%%Gaincurve: read the s,aveE to plot the gain curve
charge = 2  #nC
bc_stat = "OFF"
basedir = f'/afs/ifh.de/group/pitz/data/biaobin/sim2/2025/bc{bc_stat}_{charge}nC'
os.chdir(basedir)

# Set global styles
font = 20
plt.rcParams['font.size'] = font  # Default font size (affects labels, titles, legend)
plt.rcParams['lines.linewidth'] = 2  # Default line width
plt.rcParams['lines.markersize'] = 10  # Default marker size
plt.rcParams['legend.fontsize'] = font  # Default legend font size
plt.rcParams['xtick.labelsize'] = font  # X-axis tick label size
plt.rcParams['ytick.labelsize'] = font  # Y-axis tick label size

degl = np.arange(-16, -40 -1, -2)
# degl = [-36]

plt.figure(figsize=(20,12))
markers = ['o', 's', '^', 'v', '<', '>', 'd', 'D', 'p', '*', 'h', 'H', 'x', '+', '.', ',']
# colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'orange', 'purple', 'brown', 'pink', 'gray', 
#           '#FF5733', '#33FF57']

j=0
for deg in degl:  
    tmp = "cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{:.1f}deg/02_gen4/".format(deg)

    fname = "20seedsEnergy.txt"
    dat = np.loadtxt(tmp+fname)
    
    plt.plot(dat[::8,0],dat[::8,1],'-',marker=markers[j],label="phi2="+str(deg)+"deg")
    
    # plt.yscale('log')

    plt.legend()
    plt.xlabel('z (m)')
    plt.ylabel('<Energy> (uJ)')    
    j=j+1
plt.title(f"Q={charge} nC, BC {bc_stat}")    
plt.grid()
plt.savefig("phi2-allGainCurves-Energy"+fname+".png")


#%%%  phi2-Energy/Power with error bar

#%%%% 20seeds: get all results in loops
# get the phi2-THzEnergy txt data for each deg
profile = "flattop"

bc_statl = ["OFF","ON"]
chargel =[1,1.5,2]

seedl = np.arange(5,25,1)

if profile =="gauss":
    core=6464
    Np="5.24e+06"
    degl = np.arange(-16, -40 -1, -2)
elif profile =="flattop":
    core=6464
    Np="5.24e+06"
    degl = np.arange(0, -26-1, -2)

zenpow = "zpower" #zpower
zenpow = "zenergy"

for bc_stat in bc_statl:
    for charge in chargel:
        basedir=f'/lustre/fs22/group/pitz/biaobin/2025/{profile}_bc{bc_stat}_{charge}nC'
        os.chdir(basedir)
        
        # degl=[-16]
        for deg in degl:
            # uncooment the follwing #lines, will plot one figure for one deg with 20 shots
            plt.figure(figsize=(8,6))
            energy_20seeds = []
            for seed in seedl:  
                tmp = f"cores_016_impt_{core}0064_impz_{core}0064_Np_{Np}_phi2_{deg:.1f}deg/02_gen4/"
                
                if bc_stat == "OFF":
                    fname = f"gen4.seed{seed}.nper5.fitdeg1.out.h5"
                elif bc_stat == "ON":
                    fname = f"gen4.seed{seed}.nper10.fitdeg3.out.h5"
                # print("fname=",fname)
                pg = PostGenesis(tmp+fname, debug=0, version=4, fig_ext = '.png')
                if zenpow == "zenergy": 
                    pgzep = pg.zenergy*1e6  #uJ
                    unit = "uJ"
                    label = "THz Energy"
                elif zenpow == "zpower":
                    pgzep = pg.zpower/1e6   #MW
                    unit = "MW"
                    label = "THz power"
                
                energy_20seeds.append(pgzep)    
            
                plt.plot(pg.zplot,pgzep,"g-")
                
                plt.yscale('log')
                plt.xlabel('z (m)')
                plt.ylabel(f'{label} ({unit})')    
                
            plt.grid()
            plt.savefig(tmp+f"phi2-20seeds-GainCurves-{zenpow}.png")
            plt.close()
        
            s = pg.zplot
            aveE = np.mean(energy_20seeds,0)
            maxE = np.max(energy_20seeds,0)
            minE = np.min(energy_20seeds,0)
            sigE = np.std(energy_20seeds,0)
            
            lower_errors = np.array(aveE)-np.array(minE)
            upper_errors = np.array(maxE)-np.array(aveE)
            errors = [lower_errors, upper_errors]
            
            # plt.figure(figsize=(8,5))
            # plt.errorbar(s, aveE, yerr=errors, fmt="o-", color="gray", capsize=5, capthick=1, elinewidth=1.5)
            
            data = np.column_stack((s,aveE,minE,maxE,lower_errors,upper_errors,sigE))
            np.savetxt(tmp+f"20seeds{zenpow}.txt", data, fmt="%20.15e")
        
        #% 20seeds: (2)phi2-ave_max_E
        #---------------------------
        # run the previous section, to get the 20seeds{zenpow}.txt file first
        
        os.chdir(basedir)
        degl = np.arange(-16, -40 -1, -2)
        
        Eexit=[]
        errors_tt =[]
        sigE_exit_tt = []
        for deg in degl:
            tmp = f"cores_016_impt_{core}0064_impz_{core}0064_Np_{Np}_phi2_{deg:.1f}deg/02_gen4/"
                
            dat = np.loadtxt(tmp+f"20seeds{zenpow}.txt")
            
            aveE = dat[:,1]
            minE = dat[:,2]
            maxE = dat[:,3]
            lowErr = dat[:,4] 
            uppErr = dat[:,5]
            sigE   = dat[:,6]
            
            Eexit.append( aveE[-1] )
            errors_tt.append( [lowErr[-1],uppErr[-1]] )
            sigE_exit_tt.append( [sigE[-1],sigE[-1]])
        
        errors = np.array(errors_tt).T
        sigE_exit = np.array(sigE_exit_tt).T
            
        plt.figure(figsize=(8,6))
        # plt.plot(degl,Eexit,'o-')
        
        #p2p
        plt.errorbar(degl, Eexit, yerr=errors, fmt="o-", color="red", capsize=5, capthick=1, elinewidth=1.5)   
        #rms
        # plt.errorbar(degl, Eexit, yerr=sigE_exit, fmt="o-", color="red", capsize=5, capthick=1, elinewidth=1.5)   
        plt.grid()  
        plt.savefig(f"phi2-THz{zenpow}_bc{bc_stat}_{charge}nC.png")
        
        data = np.column_stack((degl, Eexit, errors[0], errors[1], sigE_exit[0], sigE_exit[1]))  
        np.savetxt(basedir+f"/phi2-{zenpow}Exit-Error.txt", data, fmt="%20.15e")
#%%%% 20seeds: plot together
rootdir = "/lustre/fs22/group/pitz/biaobin/2025/"

# profilel = ["gauss","flattop"]
profilel = ["gauss"]
bc_statl = ["OFF","ON"]

# bc_statl = ["OFF"]
chargel = [1,1.5,1.75,2]

zenpow = "zpower" #zpower
zenpow = "zenergy"

if zenpow == "zenergy":
    unit = "uJ"
    label="THz energy"
elif zenpow == "zpower":
    unit = "MW"
    label="THz power"

plt.figure(figsize=(8*1.5,6*1.5))
for profile in profilel:
    for bc_stat in bc_statl:
        for charge in chargel:
            
            basedir=f'/lustre/fs22/group/pitz/biaobin/2025/{profile}_bc{bc_stat}_{charge}nC'
            os.chdir(basedir)
            
            fname = f"phi2-{zenpow}Exit-Error.txt"
    
            dat = np.loadtxt(fname)
            degl = dat[:,0]
            Eexit = dat[:,1]
            
            #peak 2 peak
            lowErr = dat[:,2]
            uppErr = dat[:,3]
            errors = np.array([lowErr,uppErr])
            
            #rms, sig,sig
            # lowErr = dat[:,4]
            # uppErr = dat[:,5]
            # errors = np.array([lowErr,uppErr])
            
            plt.errorbar(degl, Eexit, yerr=errors, markersize=10 ,capsize=5, capthick=1, elinewidth=1.5,label=f'BC {bc_stat:>3}, {charge:.2f} nC')   
            plt.legend()    

plt.grid()
plt.xlabel("phi2 (deg)")
plt.ylabel(f"{label} ({unit})")
plt.savefig(f"{rootdir}{profile}-phi2-THz{zenpow}-ErrorBar-BC-{'_'.join(bc_statl)}.png")            

#%%%% seed-Energy
basedir = '/afs/ifh.de/group/pitz/data/biaobin/sim2/2025/bcON_2nC'
# basedir = '/afs/ifh.de/group/pitz/data/biaobin/sim2/2025/bcOFF_1.5nC'
# basedir = '/afs/ifh.de/group/pitz/data/biaobin/sim2/2025/bcOFF_2nC'

os.chdir(basedir)

degl = np.arange(-16, -40 -1, -2)

degl=[-28]
for deg in degl:
    # seedl = np.arange(5,26,1)
    seedl = np.arange(6,26,1)
    
    # uncooment the follwing #lines, will plot one figure for one deg with 20 shots
    plt.figure(figsize=(8,6))
    energy_20seeds = []
    for seed in seedl:  
        tmp = "cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{:.1f}deg/02_gen4/".format(deg)
        
        fname = f"gen4.seed{seed}.nper5.fitdeg1.out.h5"
        # print("fname=",fname)
        pg = PostGenesis(tmp+fname, debug=0, version=4, fig_ext = '.png')
        energy_20seeds.append(pg.zenergy[::1]*1e6)    
    
        plt.plot(seed,pg.zenergy[-1]*1e6,"^r")
        
        plt.xlabel('seed')
        plt.ylabel('Energy (uJ)')    
        
plt.grid()
plt.savefig(tmp+"/seed-Energy-2.png")
#%%% phi2-gainCurve-power
# Set global styles
font = 20
plt.rcParams['font.size'] = font  # Default font size (affects labels, titles, legend)
plt.rcParams['lines.linewidth'] = 2  # Default line width
plt.rcParams['lines.markersize'] = 10  # Default marker size
plt.rcParams['legend.fontsize'] = font  # Default legend font size
plt.rcParams['xtick.labelsize'] = font  # X-axis tick label size
plt.rcParams['ytick.labelsize'] = font  # Y-axis tick label size

os.chdir(basedir)

degl = np.arange(-16, -40 -1, -2)
# degl = [-36]

plt.figure(figsize=(20,12))
markers = ['o', 's', '^', 'v', '<', '>', 'd', 'D', 'p', '*', 'h', 'H', 'x', '+', '.', ',']
# colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'orange', 'purple', 'brown', 'pink', 'gray', 
#           '#FF5733', '#33FF57']

j=0
maxPower = []
for deg in degl:  
    tmp = "cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{:.1f}deg/02_gen4/".format(deg)

    fname = "gen4.seed5.nper5.fitdeg1.out.h5"
    pg = PostGenesis(tmp+fname, debug=0, version=4, fig_ext = '.png')
    maxPower.append(np.max(pg.zpower))

    # plt.plot(pg.zplot[::5],pg.zenergy[::5]*1e6,'-',color=colors[j],marker=markers[j],label="phi2="+str(deg)+"deg")
    plt.plot(pg.zplot[::5],pg.zpower[::5]/1e6,'-',marker=markers[j],label="phi2="+str(deg)+"deg")

    plt.yscale('log')
    
    plt.legend()
    plt.xlabel('z (m)')
    plt.ylabel('Power (MW)')    
    j=j+1
    
plt.grid()
plt.savefig("phi2-allGainCurves-Power"+fname+".png")


plt.figure(figsize=(8,6))
plt.plot(degl,np.array(maxPower)/1e6,'o-')
plt.xlabel("phi2 (deg)")
plt.ylabel("THz energy@undu exit (MW)")
plt.grid()
plt.savefig("phi2-THzPower-"+fname+".png")


#%%% phi2-Lsat
#%%%% seed5
profile = "flattop"
charge = 2  #nC
bc_stat = "ON"

basedir = f"/lustre/fs22/group/pitz/biaobin/2025/{profile}_bc{bc_stat}_{charge}nC"
os.chdir(basedir)

# degl = np.arange(-16, -40 -1, -2)
degl = np.arange(0, -26-1, -2)

# degl = [-36]

j=0
Lsat = []
for deg in degl:  
    tmp = "cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{:.1f}deg/02_gen4/".format(deg)
    
    if bc_stat == "ON":
        fname = "gen4.seed5.nper10.fitdeg3.out.h5"
    else:
        fname = "gen4.seed5.nper5.fitdeg1.out.h5"
    pg = PostGenesis(tmp+fname, debug=0, version=4, fig_ext = '.png')

    interp_func_Lsat= interp1d(pg.zenergy, pg.zplot, kind='linear', fill_value="extrapolate")
    
    maxEnergy = np.max(pg.zenergy)
    Lsat.append( interp_func_Lsat(0.9*maxEnergy) )
  
data = np.column_stack((degl,Lsat))
np.savetxt(basedir+"/seed5_degl_Lsat.txt", data, fmt="%20.15e")    

    
plt.figure(figsize=(8,6))
plt.plot(degl,Lsat,'o-')
plt.xlabel("phi2 (deg)")
plt.ylabel("sat. len. (m)")
plt.grid()
plt.savefig("phi2-Lsat.png")

#%%%% phi2-Lsat from 20seeds
# get Lsat from <aveE> gain curve

# profile = "gauss"
# bc_statl = ["OFF","ON"]
# chargel = [1,1.5,1.75,2]  #nC
# degl = np.arange(-16, -40 -1, -2)

profile = "flattop"
bc_statl = ["OFF"]
chargel = [2]  #nC
degl = np.arange(0, -26-1, -2)

rootdir = r'/lustre/fs22/group/pitz/biaobin/2025/'

for bc_stat in bc_statl:
    for charge in chargel:
        
        basedir = f"/lustre/fs22/group/pitz/biaobin/2025/{profile}_bc{bc_stat}_{charge}nC"
        os.chdir(basedir)
        
        Lsat = []
        for deg in degl:  
            tmp = "cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{:.1f}deg/02_gen4/".format(deg)
        
            fname = "20seedsEnergy.txt"
            dat = np.loadtxt(tmp+fname)
        
            interp_func_Lsat= interp1d(dat[:,1],dat[:,0], kind='linear', fill_value="extrapolate")
            
            maxEnergy = np.max(dat[:,1])
            Lsat.append( interp_func_Lsat(0.8*maxEnergy) )
            
        data = np.column_stack((degl,Lsat))
        np.savetxt(basedir+"/20seeds_ave_degl_Lsat.txt", data, fmt="%20.15e")

#%%%% 20seeds, now plot for ON and OFF
bc_statl = ["OFF","ON"]
plt.figure(figsize=(8*1.5,6*1.5))
for bc_stat in bc_statl:
    for charge in chargel:
        basedir = f"/lustre/fs22/group/pitz/biaobin/2025/{profile}_bc{bc_stat}_{charge}nC"
        fname = basedir+"/20seeds_ave_degl_Lsat.txt"
        
        dat= np.loadtxt(fname)
        
        plt.plot(dat[:,0],dat[:,1],label=f'BC {bc_stat}, {charge} nC')
        plt.legend()
    
plt.grid()
plt.xlabel("phi2 (deg)")
plt.ylabel("Sat. Length (m)") 
plt.savefig(rootdir+f"{profile}_bc{'-'.join(bc_statl)}_phi2_Lsat.png")


#%%% t-power: seed5
# didn't scan all the seeds
#%%%% (1) phi2: t-power
# profile = "gauss"
# chargel = [1, 1.5, 2]  #nC
# bc_statl = ["OFF","ON"]
# degl = np.arange(-16, -40 -1, -2)

profile = "flattop"
chargel = [1, 1.5, 2]  #nC
bc_statl = ["ON"]
degl = np.arange(0, -26-1, -2)

rootdir = r'/lustre/fs22/group/pitz/biaobin/2025/'


# Define Gaussian function
from scipy.optimize import curve_fit
def gaussian(x, a, mu, sigma):
    return a * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

key = "power"

for bc_stat in bc_statl:
    for charge in chargel:
        basedir = f"/lustre/fs22/group/pitz/biaobin/2025/{profile}_bc{bc_stat}_{charge}nC/"
        os.chdir(basedir)
        
        dat = np.loadtxt("seed5_degl_Lsat.txt")
        Lsat = dat[:,1]
        pos_sat = Lsat
        
        j=0
        sigl = []
        for deg in degl:    
            tmp = "cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{:.1f}deg/02_gen4/".format(deg)
            if bc_stat == "OFF":
                fname = "gen4.seed5.nper5.fitdeg1.out.h5"
            else:
                fname = "gen4.seed5.nper10.fitdeg3.out.h5"
                
            # pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')
            file = h5py.File(tmp+fname,'r')
            
            jj_sat = int(pos_sat[j]/3.6*file["Field"][key].shape[0])
            
            sliceNum = file["Field"][key].shape[1]
            power0 = file["Field"][key][jj_sat,:]
            t0 = np.array(list(range(0,sliceNum)))*100e-6/2.998e8 *1e12 #ps
            max_index = np.argmax(power0)
            t0 -= t0[max_index]
            
            # cut between [-20,20]
            idd = (t0>-20) & (t0<40)
            t = t0[idd]
            power = power0[idd]
            
            max_index = np.argmax(power)
            t -= t[max_index]
            
            # Initial guesses: amplitude=1, mean=0, std dev=1
            p0 = [np.max(power), 0, 1]
            # Perform curve fitting
            params, covariance = curve_fit(gaussian, t, power, p0=p0)
            power_fit = gaussian(t, *params)
            # Extract fitted parameters
            a_fit, mu_fit, sigma_fit = params
            print(f"Fitted Parameters: Amplitude={a_fit:.3f}, Mean={mu_fit:.3f}, Sigma={sigma_fit:.3f}")
            sigl.append(sigma_fit)
            
            plt.figure(figsize=(8,6))
            plt.plot(t[::1],power[::1],label="phi2="+str(deg)+"deg")
            plt.plot(t, power_fit,'--',label='fit, $\sigma_t$=%.2f ps' % (sigma_fit))
            
            plt.xlabel('t (ps)')
            plt.ylabel(key+"@sat (W)")
            plt.legend()
            plt.grid()
            plt.savefig(basedir+f"/t_power_fitting_phi2={deg}deg")
            plt.close()
            
            j=j+1
        
        dat = np.column_stack((degl,sigl))
        np.savetxt(basedir+"degl_sigl.txt",dat)

#%%%% (2) plot together

def plot_phi2_rmsPul(bc_stat):
    plt.figure(figsize=(8,6))
    for charge in chargel:
        basedir = f'{rootdir}{profile}_bc{bc_stat}_{charge}nC'
        fname = basedir+"/degl_sigl.txt"
        
        dat= np.loadtxt(fname)
    
        plt.plot(dat[:,0],dat[:,1],label=f'BC {bc_stat}, {charge} nC')
        
    plt.grid()
    plt.legend()
    plt.xlabel("phi2 (deg)")
    plt.ylabel("rms THz pulse Len. (ps)")  
    plt.ylim([1,None])
      
    plt.savefig(rootdir+f"{profile}_phi2_PulseLen_bc{bc_stat}.png")


plot_phi2_rmsPul("ON")
# plot_phi2_rmsPul("OFF")



#%%% t-power: scan 20 seeds for each phi2
rootdir = r'/lustre/fs22/group/pitz/biaobin/2025/'

profile = "flattop"
chargel = [1, 1.5, 2]  #nC
bc_statl = ["OFF","ON"]
# degl = np.arange(-16, -40 -1, -2)
degl = np.arange(0, -26 -1, -2)
seedl = np.arange(5,25,1)

# profile = "gauss"
# chargel = [2]  #nC
# bc_statl = ["OFF"]
# degl = [-32]
# seedl = [5]

# Define Gaussian function
from scipy.optimize import curve_fit
def gaussian(x, a, mu, sigma):
    return a * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

for bc_stat in bc_statl:
    for charge in chargel:
    
        basedir = f'/lustre/fs22/group/pitz/biaobin/2025/{profile}_bc{bc_stat}_{charge}nC/'
        os.chdir(basedir)
        
        dat = np.loadtxt("20seeds_ave_degl_Lsat.txt")
        Lsat = dat[:,1]
        pos_sat = dat[:,1] 
        
        j=0
        key = "power"
        
        min_sigt=[]
        max_sigt=[]
        ave_sigt=[]
        for deg in degl:    
            tmp = "cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{:.1f}deg/02_gen4/".format(deg)
            
            plt.figure(figsize=(20/2,12/2))            
            sigl = []
            for seed in seedl:
                if bc_stat == "OFF":
                    fname = f"gen4.seed{seed}.nper5.fitdeg1.out.h5"
                else:
                    fname = f"gen4.seed{seed}.nper10.fitdeg3.out.h5"
                
                # pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')
                file = h5py.File(tmp+fname,'r')
                
                jj_sat = int(pos_sat[j]/3.6*file["Field"][key].shape[0])
                
                sliceNum = file["Field"][key].shape[1]
                power0 = file["Field"][key][jj_sat,:]
                t0 = np.array(list(range(0,sliceNum)))*100e-6/2.998e8 *1e12 #ps
                max_index = np.argmax(power0)
                t0 -= t0[max_index]
                
                # cut between [-20,20]
                # idd = (t0>-100) & (t0<100)
                idd = (t0>-20) & (t0<40)
                t = t0[idd]
                power = power0[idd]
                
                max_index = np.argmax(power)
                t -= t[max_index]
                
                # Initial guesses: amplitude=1, mean=0, std dev=1
                p0 = [np.max(power), 0, 1]
                # Perform curve fitting
                params, covariance = curve_fit(gaussian, t, power, p0=p0)
                power_fit = gaussian(t, *params)
                # Extract fitted parameters
                a_fit, mu_fit, sigma_fit = params
                print(f"Fitted Parameters: Amplitude={a_fit:.3f}, Mean={mu_fit:.3f}, Sigma={sigma_fit:.3f}")
                sigl.append(sigma_fit)
                
                
                plt.plot(t[::1],power[::1],'-',marker=markers[j])
                plt.plot(t, power_fit,'--',label='fit, $\sigma_t$=%.2f ps' % (sigma_fit))
            plt.title("phi2="+str(deg)+"deg")   
            plt.xlabel('t (ps)')
            plt.ylabel(key+"@sat (W)")
            plt.legend()
            plt.grid()
            plt.title(f"{profile}, BC {bc_stat}, {charge:.2f} nC @phi2={deg} deg")
            plt.savefig(basedir+f"/bc_{bc_stat}_t_power_fitting_20seeds_phi2={deg}deg")
            # plt.close()
            
            min_sigt.append(np.min(sigl))
            max_sigt.append(np.max(sigl))
            ave_sigt.append(np.mean(sigl))
            
            j=j+1
            
        
        dat = np.column_stack((degl,ave_sigt,min_sigt,max_sigt))
        np.savetxt(basedir+"/20seeds_degl_ave_min_max_sigt.txt",dat)

#%%%% plot together

bc_statl = ["OFF","ON"]
# bc_statl = ["OFF"]

#---------------------------------------  
plt.figure(figsize=(8,6))  

for bc_stat in bc_statl:
    for charge in chargel:
        basedir = rootdir+f'{profile}_bc{bc_stat}_{charge}nC/'
        fname = "20seeds_degl_ave_min_max_sigt.txt"
        
        dat = np.loadtxt(basedir+fname)
        degl  = dat[:,0]
        avePL = dat[:,1]
        minPL = dat[:,2]
        maxPL = dat[:,3]
        lowErr= np.array(avePL)-np.array(minPL)
        uppErr= np.array(maxPL)-np.array(avePL)
    
        errors = np.array([lowErr,uppErr])
        
        plt.errorbar(degl, avePL, yerr=errors,markersize=20, \
                      capsize=5, capthick=2, elinewidth=1.5,label=f'BC {bc_stat:<3}, {charge:.2f} nC')   
        plt.legend()    
    
plt.xlabel("phi2 (deg)")
plt.ylabel("rms THz pulse Len. (ps)")
plt.grid()
plt.savefig(rootdir+f"{profile}_phi2-PulseLen-20seeds-BC{'-'.join(bc_statl)}.png")


#%%% phi2: t-intensityFarField
font = 20
plt.rcParams['font.size'] = font  # Default font size (affects labels, titles, legend)
plt.rcParams['lines.linewidth'] = 2  # Default line width
plt.rcParams['lines.markersize'] = 10  # Default marker size
plt.rcParams['legend.fontsize'] = font  # Default legend font size
plt.rcParams['xtick.labelsize'] = font  # X-axis tick label size
plt.rcParams['ytick.labelsize'] = font  # Y-axis tick label size

basedir = '/afs/ifh.de/group/pitz/data/biaobin/sim2/2025/bcOFF_1nC'
# basedir = '/afs/ifh.de/group/pitz/data/biaobin/sim2/2025/bcOFF_1.5nC'
# basedir = '/afs/ifh.de/group/pitz/data/biaobin/sim2/2025/bcOFF_2nC'
os.chdir(basedir)

degl = np.arange(-16, -40 -1, -2)
# degl = [-36]

plt.figure(figsize=(20,12))
markers = ['o', 's', '^', 'v', '<', '>', 'd', 'D', 'p', '*', 'h', 'H', 'x', '+', '.', ',']

dat = np.loadtxt("20seeds_ave_degl_Lsat.txt")
Lsat = dat[:,1]
pos_sat = dat[:,1] 
j=0

key = 'intensity-farfield'
for deg in degl:    
    tmp = "cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{:.1f}deg/02_gen4/".format(deg)
    fname = "gen4.seed5.nper5.fitdeg1.out.h5"
    
    # pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')
    file = h5py.File(tmp+fname,'r')
    
    jj_sat = int(pos_sat[j]/3.6*file["Field"][key].shape[0])
    
    sliceNum = file["Field"][key].shape[1]
    power = file["Field"][key][jj_sat,:]
    
    max_index = np.argmax(power)
    t = np.array(list(range(0,sliceNum)))*100e-6/2.998e8 *1e12 #ps
    t -= t[max_index]
    
    plt.plot(t[::1],power[::1],'-',marker=markers[j],label="phi2="+str(deg)+"deg")
    
    plt.xlabel('t (ps)')
    plt.ylabel(key+"@sat (arb.)")
    plt.legend()
    j=j+1
plt.grid()
# plt.xlim([80,120])
plt.savefig("phi2-t-"+key+"-"+fname+".png")



#%%% phi2: t-nearFarField
font = 20
plt.rcParams['font.size'] = font  # Default font size (affects labels, titles, legend)
plt.rcParams['lines.linewidth'] = 2  # Default line width
plt.rcParams['lines.markersize'] = 10  # Default marker size
plt.rcParams['legend.fontsize'] = font  # Default legend font size
plt.rcParams['xtick.labelsize'] = font  # X-axis tick label size
plt.rcParams['ytick.labelsize'] = font  # Y-axis tick label size

os.chdir(basedir)

degl = np.arange(-16, -40 -1, -2)
# degl = [-36]

plt.figure(figsize=(20,12))
markers = ['o', 's', '^', 'v', '<', '>', 'd', 'D', 'p', '*', 'h', 'H', 'x', '+', '.', ',']

pos_sat = Lsat   #get Lsat first

j=0

key = 'intensity-nearfield'
for deg in degl:    
    tmp = "cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{:.1f}deg/02_gen4/".format(deg)
    fname = "gen4.seed5.nper5.fitdeg1.out.h5"
    
    # pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')
    file = h5py.File(tmp+fname,'r')
    
    jj_sat = int(pos_sat[j]/3.6*file["Field"][key].shape[0])
    
    sliceNum = file["Field"][key].shape[1]
    power = file["Field"][key][jj_sat,:]
    
    max_index = np.argmax(power)
    t = np.array(list(range(0,sliceNum)))*100e-6/2.998e8 *1e12 #ps
    t -= t[max_index]
    
    plt.plot(t[::1],power[::1],'-',marker=markers[j],label="phi2="+str(deg)+"deg")
    
    plt.xlabel('t (ps)')
    plt.ylabel(key+"@sat (arb.)")
    plt.legend()
    j=j+1
plt.grid()
plt.savefig("phi2-t-"+key+"-"+fname+".png")

#%%% phi2-spectrum
rootdir = f"/lustre/fs22/group/pitz/biaobin/2025/"
os.chdir(rootdir)

profile = "flattop"
# bc_stat = "ON"
bc_stat = "OFF"

chargel = [1,1.5,2]


# degl = [-32]

plt.figure(figsize=(8,6))
markers = ['o', 's', '^', 'v', '<', '>', 'd', 'D', 'p', '*', 'h', 'H', 'x', '+', '.', ',']

deg=-14
for charge in chargel:
    basedir = f"/lustre/fs22/group/pitz/biaobin/2025/{profile}_bc{bc_stat}_{charge}nC/"
    os.chdir(basedir)
    
    # dat = np.loadtxt("20seeds_ave_degl_Lsat.txt")
    dat = np.loadtxt("seed5_degl_Lsat.txt")
    degll = dat[:,0]
    pos_satl = dat[:,1] 
    
    tmp = "cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_{:.1f}deg/02_gen4/".format(deg)
    
    if bc_stat =="OFF":
        fname = basedir+tmp+"gen4.seed5.nper5.fitdeg1.out.h5"
    else:
        fname = basedir+tmp+"gen4.seed5.nper10.fitdeg3.out.h5"
    
    pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')
   
    interp_func_posSat= interp1d(degll, pos_satl, kind='linear', fill_value="extrapolate")
    pos_sat = interp_func_posSat(deg)
    print(f"pos sat@{deg}deg=",pos_sat)
    
    jj_sat = int(pos_sat/3.6*pg.zplot.shape[0])
    wavelen, amp = pg.get_spectrum(pg.zplot[jj_sat])
    
    # plt.plot(wavelen*1e6,amp/np.max(amp),'-',label=f"{charge} nC")
    plt.plot(wavelen*1e6,amp,'-',label=f"{charge} nC")
    plt.legend()
    plt.xlabel('wavelength (um)')
    plt.ylabel('intensity')
    
plt.grid()
plt.title(f"{profile}, BC {bc_stat:<3}, phi2={deg} deg")
plt.xlim([60,150])

plt.savefig(f"{rootdir}{profile}_BC{bc_stat}_phi2-THZ-spectrum@sat.png")

#%%%% single file spectrum
font = 20
plt.rcParams['font.size'] = font  # Default font size (affects labels, titles, legend)
plt.rcParams['lines.linewidth'] = 2  # Default line width
plt.rcParams['lines.markersize'] = 10  # Default marker size
plt.rcParams['legend.fontsize'] = font  # Default legend font size
plt.rcParams['xtick.labelsize'] = font  # X-axis tick label size
plt.rcParams['ytick.labelsize'] = font  # Y-axis tick label size

#------------------------------------------
basedir=r'/afs/ifh.de/group/pitz/data/biaobin/sim1/2025/00_sase_chicane/Q_1_1.5_2nC/bcOFF_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/'
os.chdir(basedir)
fname = "gen4.seed5.nper5.fitdeg1.fixP0.out.h5"
pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

pos_sat = 3.0
jj_sat = int(pos_sat/3.6*pg.zplot.shape[0])
wavelen, amp = pg.get_spectrum(pg.zplot[jj_sat])

plt.plot(wavelen*1e6,amp,'--',label="phi2="+str(deg)+"deg")
plt.legend()
plt.xlabel('wavelength (um)')
plt.ylabel('intensity') 
plt.show()

#------------------------------------------
# plt.figure(figsize=(8,5))
# deg = -32
# basedir = '/afs/ifh.de/group/pitz/data/biaobin/sim1/2024/SASE_Chicane/20241202-NonCompressor-ScanPhase/cores_016_impt_32320064_impz_32320064_Np_1.00e+07_phi2_-32.0deg/02_gen4/'
# os.chdir(basedir)
# fname = "gen4.out.h5"
# pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

# pos_sat = 3.0
# jj_sat = int(pos_sat/3.6*pg.zplot.shape[0])
# wavelen, amp = pg.get_spectrum(pg.zplot[jj_sat])

# plt.plot(wavelen*1e6,amp,'-',label="phi2="+str(deg)+"deg")
# plt.legend()
# plt.xlabel('wavelength (um)')
# plt.ylabel('intensity') 
# plt.show()

#%%% given phi2, animation


# Set global styles
font = 12
plt.rcParams['font.size'] = font  # Default font size (affects labels, titles, legend)
plt.rcParams['lines.linewidth'] = 2  # Default line width
plt.rcParams['lines.markersize'] = 10  # Default marker size
plt.rcParams['legend.fontsize'] = font  # Default legend font size
plt.rcParams['xtick.labelsize'] = font  # X-axis tick label size
plt.rcParams['ytick.labelsize'] = font  # Y-axis tick label size

# workdir ='/lustre/fs22/group/pitz/biaobin/2025/test_flattop_rs_3ps_22psTop/flattop_bcOFF_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-8.0deg/02_gen4/copenChirp'
# os.chdir(workdir)
# fname = "../gen4.seed5.nper5.fitdeg1.out.h5"

# test for 170MeV beam
# path = '/lustre/fs22/group/pitz/biaobin/2025/gauss_bcON_1nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-32.0deg/02_gen4/test_170MeV'
# os.chdir(path)
# fname = "gen4.seed5.nper10.fitdeg3.out.h5"

fname = "gen4.out.h5"

pg = PostGenesis(fname, debug=0, version=4, fig_ext = '.png')

fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(nrows = 5, figsize = (4, 10))
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

step = 10
for i in np.arange(0, ppower.shape[0], step):
# for i in [20]:
    
    # ax1.plot(p_mid[i]/p_mid[i].max(), 'r-')
    # ax12.plot(ppower[i]/ppower[i].max(), 'b-')
    ax1.plot(p_mid[i]/p_mid[i].max(), 'r-',label='intensity-farfield')
    ax1.plot(ppower[i]/ppower[i].max(), 'b-',label='power')
    
    # title = str.format('z = %.3f m,%0d undu.' % (pg.zplot[i], (pg.zplot[i]-0.105)/0.03))
    title = str.format('z = %.3f m' % (pg.zplot[i]))

    ax1.set_title(title)
    ax1.set_xlabel(r'# of slice')
    ax1.set_ylabel(r'On-axis power (au)')
    #ax12.set_ylabel(r'Power (au)')
    #ax.set_ylim(0, 2500)
    # ax1.grid()
    ax1.set_xlim(0, 125)
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
    # ax2.set_ylim(25, 40)
    # ax2.set_xlim(0, 125)
    ax2.set_xlim(0, 125)
    ax2.legend()
    
    DE = (gamma[i]-gamma[0])*pg.current
    #ax3.plot(DE/np.max(np.abs(DE)), 'r-')
    ax3.plot(DE, 'r-')
    ax3.set_xlim(0, 125)
    # ax3.set_ylim(-2, 2)
    ax3.set_xlabel('# of slice')
    ax3.set_ylabel(r'Energy loss')
    
    ax32.plot(pg.current, 'k-',label='current prof.')
    #ax32.set_ylabel(r'Current (A)')
    ax32.legend()
    
    ax4.plot(bunching[i], 'r-')
    ax4.set_xlabel('# of slice')
    ax4.set_ylabel('Bunching')
    # ax4.set_ylim(1e-5, 1)
    ax4.set_yscale('log')
    ax4.set_xlim(0, 125)
    
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
    fig.savefig('00_z_%.2fm' % pg.zplot[i] + '.png')    
    
    # plt.pause(1)
    if i<ppower.shape[0]-step:
        ax1.cla()
        ax12.cla()
        ax2.cla()
        ax3.cla()
        ax32.cla()
        ax4.cla()
        ax5.cla()

# os.chdir(basedir)

#%%% initial beam phase space
workdir ='/lustre/fs22/group/pitz/biaobin/2025/test_flattop_rs_3ps_22psTop/flattop_bcOFF_2nC/cores_016_impt_64640064_impz_64640064_Np_5.24e+06_phi2_-8.0deg/02_gen4/copenChirp'
os.chdir(workdir)

fname = "../scan.seed5.nper5.fitdeg1.out.par.h5"

f = h5py.File(fname,'r')

dlamda = 0 
lamdas = 1 
plt.figure(figsize=(8,6))
for j in range(f['slicecount'][0]):
# for j in np.arange(10,15):
    print("j=",j)
    if j+1 < 10: 
        name = 'slice000'+'00'+str(j+1)
    elif j+1 < 100:
        name = 'slice000'+'0'+str(j+1)
    elif j+1 < 1000:
        name = 'slice000'+str(j+1)
    plt.plot(f[name]['theta'][:]/(2*np.pi)*lamdas+dlamda,f[name]['x'][:],'.')
    dlamda = dlamda +lamdas
plt.xlabel('z/lambda0')
plt.ylabel('gama')

#%%%% astra dist
dat = np.loadtxt("fort.1032.001")

z = dat[:,2]*1e3 #mm
dgam = dat[:,5] 






