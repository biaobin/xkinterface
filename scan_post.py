import h5py
import xkinterface.startup as stp
import os
import matplotlib.pyplot as plt
import numpy as np
from xkinterface.interface import postgen4

stp.setplot()

try:
    hid.close()
except:
    pass

#%% test, single file

# path = '/lustre/fs25/group/pitz/biaobin/202510_idealMachine/parabolic_slscON_llscON_2/Ipeak190A_Lbuncht9.70ps/02_gen4'

# fname = "/gen4.out.h5"
# hid = h5py.File(path+fname,'r')

# current = postgen4.plot_BeamCurrentProfile(hid)


#%% now imshow for Lbunchtl-Ipeak VS pulse energy
unduu = "Planar"
profile="parabolic"
# profile="flattop"

slsc="ON"
llsc="ON"

rootdir = f"/lustre/fs25/group/pitz/biaobin/202511_idealMachine/{unduu}_{profile}_slsc{slsc}_llsc{llsc}"
os.chdir(rootdir)

# Lbunchtl = np.arange(0.1e-12, 30e-12, 1e-12)  # [s]
# Ipeakl= np.arange(60,600,20)  #A

# scan a smaller range
# Lbunchtl = np.arange(0.1e-12, 20e-12, 1e-12)  # [s]
# Ipeakl= np.arange(60,600,20)  #A

# scan a smaller range
Lbunchtl = np.arange(0.1e-12, 15e-12, 1e-12)  # [s]
Ipeakl= np.arange(60,400,20)  #A

energyl = []
peakPowerl=[]
chargel=[]
for Lbuncht in Lbunchtl:
    for Ipeak in Ipeakl:

        fname = f"Ipeak{Ipeak:.0f}A_Lbuncht{Lbuncht*1e12:.2f}ps/02_gen4/gen4.out.h5"
        
        print("fname=",fname)
        with h5py.File(fname,'r') as hid:
            tmp1 = hid['Field']['Global']['energy'][-1]
            tmp2 = np.max(hid["Field"]["power"][:,:])
            tmp3=np.sum(hid["Beam"]["current"][0,:])*hid["Global"]["lambdaref"][0]/2.998e8 *1e9 #nC
            
        energyl.append(tmp1)
        peakPowerl.append(tmp2)
        chargel.append(tmp3)

#%%% dump the results        
import json
data = {
    "Lbunchtl": Lbunchtl.tolist(),
    "Ipeakl": Ipeakl.tolist(),
    "energyl": [float(x) for x in energyl],
    "peakPowerl": [float(x) for x in peakPowerl],
    "chargel": [float(x) for x in chargel],
}

with open("fel_results.json", "w") as f:
    json.dump(data, f, indent=2)        
        
# load it again
with open("fel_results.json", "r") as f:
    data = json.load(f)

Lbunchtl = np.array(data["Lbunchtl"])
Ipeakl = np.array(data["Ipeakl"])
energyl = np.array(data["energyl"])
peakPowerl = np.array(data["peakPowerl"])
chargel = np.array(data["chargel"])

        
#%%% get the Lbunch=10.1 ps, Q=0.5, 1.0, 1.5, ...  => peakPower & energy
# from scipy.interpolate import RegularGridInterpolator

# Lbuncht = Lbunchtl[-1] 
# c = 2.998e8

# tmpQ = np.arange(0.5, np.max(chargel),0.5)*1e-9
# if profile=='parabolic':
#     sigz = Lbuncht/6*c

#     Ipeak2 = tmpQ/(13.3422*sigz*1e-9)
# elif profile=='flattop':
#     Ipeak2 = tmpQ/(0.9*Lbuncht)
    
# print(f"Ipeak={Ipeak2} A")  

# # find the value at (Lbuncht, Ipeak) 
# energyl2 = np.array(energyl).reshape(len(Lbunchtl), len(Ipeakl))
# peakPowerl2 = np.array(peakPowerl).reshape(len(Lbunchtl), len(Ipeakl))
# chargel2 = np.array(chargel).reshape(len(Lbunchtl), len(Ipeakl))

# energy_interp = RegularGridInterpolator((Lbunchtl, Ipeakl), energyl2)
# power_interp  = RegularGridInterpolator((Lbunchtl, Ipeakl), peakPowerl2)
# charge_interp = RegularGridInterpolator((Lbunchtl, Ipeakl), chargel2)


# Lb_target = Lbuncht

# E_interp = []
# P_interp = []
# Q_interp = []
# for Ip_target in Ipeak2:
#     tmp = energy_interp((Lb_target, Ip_target))
    # E_interp.append(tmp.item())
    # P_interp.append(power_interp((Lb_target, Ip_target)).item())
    # Q_interp.append(charge_interp((Lb_target, Ip_target)).item())            
        
# print(f"Interpolated at Lbuncht={Lb_target*1e12:.2f} ps, Ipeak={Ip_target} A:")
# print(f"  Energy = {E_interp*1e6 :.2f} uJ")
# print(f"  PeakPower = {P_interp/1e6 :.2f} MW")
# print(f"  Charge = {Q_interp} nC")

#%%%%
from scipy.interpolate import RegularGridInterpolator

Lbuncht = 10.1e-12 
c = 2.998e8

tmpQ = 2e-9
if profile=='parabolic':
    sigz = Lbuncht/6*c

    Ipeak = tmpQ/(13.3422*sigz*1e-9)
elif profile=='flattop':
    Ipeak = tmpQ/(0.9*Lbuncht)
    
print(f"Ipeak={Ipeak} A")  

# find the value at (Lbuncht, Ipeak) 
energyl2 = np.array(energyl).reshape(len(Lbunchtl), len(Ipeakl))
peakPowerl2 = np.array(peakPowerl).reshape(len(Lbunchtl), len(Ipeakl))
chargel2 = np.array(chargel).reshape(len(Lbunchtl), len(Ipeakl))

energy_interp = RegularGridInterpolator((Lbunchtl, Ipeakl), energyl2)
power_interp  = RegularGridInterpolator((Lbunchtl, Ipeakl), peakPowerl2)
charge_interp = RegularGridInterpolator((Lbunchtl, Ipeakl), chargel2)

Lb_target = Lbuncht
Ip_target = Ipeak

E_interp = energy_interp((Lb_target, Ip_target))
P_interp = power_interp((Lb_target, Ip_target))
Q_interp = charge_interp((Lb_target, Ip_target))              
        
print(f"Interpolated at Lbuncht={Lb_target*1e12:.2f} ps, Ipeak={Ip_target} A:")
print(f"  Energy = {E_interp*1e6 :.2f} uJ")
print(f"  PeakPower = {P_interp/1e6 :.2f} MW")
print(f"  Charge = {Q_interp} nC")





#%%% energy                 
energyl = np.array(energyl)
n_Lbuncht = len(Lbunchtl)
n_I = len(Ipeakl)  

Z = energyl.reshape(n_Lbuncht, n_I).T*1e3 #mJ

# prepare axes in physical units
Lbuncht_ps = Lbunchtl * 1e12    # in ps
I_A = Ipeakl             # in A

# image extent: [x_min, x_max, y_min, y_max]
extent = [Lbuncht_ps[0], Lbuncht_ps[-1], I_A[0], I_A[-1]]

plt.figure()
im = plt.imshow(Z, origin='lower', aspect='auto', extent=extent)
plt.colorbar(im, label='pulse energy (mJ)')
plt.xlabel(r'bunch length (ps)')
plt.ylabel('Ipeak (A)')
# plt.title('Energy vs Ipeak and sigma_t')
plt.tight_layout()

stp.savefig("scan_energy")

#%%% 2D-energy, with contour lines
energyl = np.array(energyl)
chargel = np.array(chargel)
n_Lbuncht = len(Lbunchtl)
n_I = len(Ipeakl)  

Z = energyl.reshape(n_Lbuncht, n_I).T*1e3 #mJ
Zc = chargel.reshape(n_Lbuncht, n_I).T            # nC

# prepare axes in physical units
Lbuncht_ps = Lbunchtl * 1e12  # in ps
I_A = Ipeakl             # in A

# image extent: [x_min, x_max, y_min, y_max]
extent = [Lbuncht_ps[0], Lbuncht_ps[-1], I_A[0], I_A[-1]]

plt.figure()
im = plt.imshow(Z, origin='lower', aspect='auto', extent=extent)
plt.colorbar(im, label='pulse energy (mJ)')
plt.xlabel(r'bunch length (ps)')
plt.ylabel('Ipeak (A)')
# plt.title('Energy vs Ipeak and sigma_t')

# # contour for total charge
# # ---------------------------
# X, Y = np.meshgrid(Lbuncht_ps, I_A)
# levels = np.linspace(np.nanmin(Zc), np.nanmax(Zc), 14)
# cs = plt.contour(X, Y, Zc, levels=levels, colors='white', linewidths=2.0)
# plt.clabel(cs, inline=True, fontsize=15, fmt='%.1f nC')

# contour for total energy
# ---------------------------
X, Y = np.meshgrid(Lbuncht_ps, I_A)

levels = np.linspace(np.nanmin(Z), np.nanmax(Z), 30)
# vmin = np.nanmin(Z)
# vmax = np.nanmax(Z)
# def smooth_scale(x):
#     return np.tanh(2*x - 1) / np.tanh(1) / 2 + 0.5
# t = np.linspace(0, 1, 10)
# levels = vmin + (vmax - vmin) * smooth_scale(t)

cs = plt.contour(X, Y, Z, levels=levels, colors='gray', linewidths=1.5,linestyles='--')
plt.clabel(cs, inline=True, fontsize=16, fmt='%.2f')
# plt.clabel(cs, inline=True, fontsize=15, fmt='%.1f mJ',manual=True)


plt.tight_layout()

#%%% 2D-energy, with contour plot for charge and energy
energyl = np.array(energyl)
chargel = np.array(chargel)
n_Lbuncht = len(Lbunchtl)
n_I = len(Ipeakl)  

Z = energyl.reshape(n_Lbuncht, n_I).T*1e3 #mJ
Zc = chargel.reshape(n_Lbuncht, n_I).T            # nC

# prepare axes in physical units
Lbuncht_ps = Lbunchtl * 1e12  # in ps
I_A = Ipeakl             # in A

# image extent: [x_min, x_max, y_min, y_max]
extent = [Lbuncht_ps[0], Lbuncht_ps[-1], I_A[0], I_A[-1]]

plt.figure()
im = plt.imshow(Z, origin='lower', aspect='auto', extent=extent)
plt.colorbar(im, label='pulse energy (mJ)')
plt.xlabel(r'bunch length (ps)')
plt.ylabel('Ipeak (A)')
# plt.title('Energy vs Ipeak and sigma_t')

# contour for total charge
# ---------------------------
X, Y = np.meshgrid(Lbuncht_ps, I_A)
# levels = np.linspace(np.nanmin(Zc), np.nanmax(Zc), 14)
# levels=np.arange(0.5, 6.5, 0.5)
levels=np.arange(0.5,int(np.nanmax(Zc)),0.5)

cs = plt.contour(X, Y, Zc, levels=levels, colors='gray', linewidths=1.5,linestyles='--')

plt.clabel(cs, inline=True, fontsize=15, fmt='%.1f nC')

# contour for total energy
# ---------------------------
X, Y = np.meshgrid(Lbuncht_ps, I_A)

levels = np.linspace(np.nanmin(Z), np.nanmax(Z), 15)
# levels1 = np.linspace(np.nanmin(Z), 0.3*np.nanmax(Z), 5)
# levels2 = np.linspace(0.3*np.nanmax(Z),np.nanmax(Z), 5)
# levels = np.append(levels1,levels2[1:])

# vmin = np.nanmin(Z)
# vmax = np.nanmax(Z)
# def smooth_scale(x):
#     return np.tanh(2*x - 1) / np.tanh(1) / 2 + 0.5
# t = np.linspace(0, 1, 15)
# levels = vmin + (vmax - vmin) * smooth_scale(t)

cs = plt.contour(X, Y, Z, levels=levels, colors='white', linewidths=1.5,linestyles='-')
plt.clabel(cs, inline=True, fontsize=15, fmt='%.2f')
# plt.clabel(cs, inline=True, fontsize=15, fmt='%.1f mJ',manual=True)


plt.tight_layout()

stp.savefig("scan_energy")


#%%% peak power                
peakPowerl = np.array(peakPowerl)
n_Lbuncht = len(Lbunchtl)
n_I = len(Ipeakl)  

Z = peakPowerl.reshape(n_Lbuncht, n_I).T/1e6 #MW

# prepare axes in physical units
Lbuncht_ps = Lbunchtl * 1e12   # in ps
I_A = Ipeakl             # in A

# image extent: [x_min, x_max, y_min, y_max]
extent = [Lbuncht_ps[0], Lbuncht_ps[-1], I_A[0], I_A[-1]]

plt.figure()
im = plt.imshow(Z, origin='lower', aspect='auto', extent=extent)
plt.colorbar(im, label='peak power (MW)')
plt.xlabel(r'bunch length (ps)')
plt.ylabel('Ipeak (A)')
# plt.title('Energy vs Ipeak and sigma_t')
plt.tight_layout()

stp.savefig("scan_peakPower")

#%%% charge 
chargel = np.array(chargel)
n_Lbuncht = len(Lbunchtl)
n_I = len(Ipeakl)  

Z = chargel.reshape(n_Lbuncht, n_I).T #nC

# prepare axes in physical units
Lbuncht_ps = Lbunchtl * 1e12  # in ps
I_A = Ipeakl             # in A

# image extent: [x_min, x_max, y_min, y_max]
extent = [Lbuncht_ps[0], Lbuncht_ps[-1], I_A[0], I_A[-1]]

plt.figure()
im = plt.imshow(Z, origin='lower', aspect='auto', extent=extent)
plt.colorbar(im, label='total charge (nC)')
plt.xlabel(r'bunch length (ps)')
plt.ylabel('Ipeak (A)')
# plt.title('Energy vs Ipeak and sigma_t')
plt.tight_layout()

stp.savefig("scan_charge")


#%%% 2D-peak_power, with contour plot for charge and energy
import numpy as np

peakPowerl = np.array(peakPowerl)
chargel = np.array(chargel)
n_Lbuncht = len(Lbunchtl)
n_I = len(Ipeakl)  

Z = peakPowerl.reshape(n_Lbuncht, n_I).T/1e6 #MW
Zc = chargel.reshape(n_Lbuncht, n_I).T            # nC

# prepare axes in physical units
Lbuncht_ps = Lbunchtl * 1e12  # in ps
I_A = Ipeakl             # in A

# image extent: [x_min, x_max, y_min, y_max]
extent = [Lbuncht_ps[0], Lbuncht_ps[-1], I_A[0], I_A[-1]]

plt.figure()
im = plt.imshow(Z, origin='lower', aspect='auto', extent=extent)
plt.colorbar(im, label='peak power (MW)')
plt.xlabel(r'bunch length (ps)')
plt.ylabel('Ipeak (A)')
# plt.title('Energy vs Ipeak and sigma_t')

# contour for total charge
# ---------------------------
X, Y = np.meshgrid(Lbuncht_ps, I_A)
# levels = np.linspace(np.nanmin(Zc), np.nanmax(Zc), 15)
# levels=np.arange(0.5, 6.5, 0.5)
levels=np.arange(0.5,int(np.nanmax(Zc)),0.5)

cs = plt.contour(X, Y, Zc, levels=levels, colors='gray', linewidths=1.5,linestyles='--')

plt.clabel(cs, inline=True, fontsize=15, fmt='%.1f nC')

# contour for peak power
# ---------------------------
X, Y = np.meshgrid(Lbuncht_ps, I_A)
levels = np.linspace(np.nanmin(Z), np.nanmax(Z), 12)


# vmin = np.nanmin(Z)
# vmax = np.nanmax(Z)
# def smooth_scale(x):
#     return np.tanh(2*x - 1) / np.tanh(1) / 2 + 0.5
# t = np.linspace(0, 1, 10)
# levels = vmin + (vmax - vmin) * smooth_scale(t)

cs = plt.contour(X, Y, Z, levels=levels, colors='white', linewidths=1.5,linestyles='-')
plt.clabel(cs, inline=True, fontsize=15, fmt='%.1f')
# plt.clabel(cs, inline=True, fontsize=15, fmt='%.1f mJ',manual=True)


plt.tight_layout()

stp.savefig("scan_2D_peakPower")





