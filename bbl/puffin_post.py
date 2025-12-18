import h5py
import xkinterface.startup as stp
import os
import matplotlib.pyplot as plt
import numpy as np
import tables
import xkinterface.startup as stp
stp.setplot()


stp.setplot(figsize=(12,8))


# path= "/mnt/f/simu_2025/202508_Puffin/PhyOfPlasmasV19pp093119/fig7a"
# path= "/mnt/f/simu_2025/202508_Puffin/PhyOfPlasmasV19pp093119/fig7"
# path = stp.winpath(r"F:\simu_2025\202508_Puffin\test4_fromxiangkun")

# os.chdir(path)

fname = "lcls_integrated_70.h5"

# fname = "lcls_integrated_101.h5"
hid = h5py.File(fname)

h5f = tables.open_file(fname, mode='r')
zD = h5f.root.runInfo._v_attrs.zTotal

#%% power profile in SI unit
# nz2 =  h5f.root.runInfo._v_attrs.nZ2
# dz2 =  h5f.root.runInfo._v_attrs.sLengthOfElmZ2
# lc  =  h5f.root.runInfo._v_attrs.Lc

spos = hid["runInfo"].attrs["zTotal"]

nz2 = hid["runInfo"].attrs["nZ2"]
dz2 = hid["runInfo"].attrs["sLengthOfElmZ2"]
lc  = hid["runInfo"].attrs["Lc"]

c0  =  2.99792458e8

lenz2 = (nz2-1) *dz2
z2 = (np.arange(0,nz2)) *dz2
s = z2 *lc
t = s / c0

power = hid["powerSI"][:]

# plt.figure()
# plt.title(f"s={spos:.2f} (m)")
# plt.plot(t,power,'-')
# plt.xlabel('t (s)')
# plt.ylabel("power (w)")

# plt.figure()

fig, ax1 = plt.subplots()

plt.title(f"s={spos:.2f} (m)")
ax1.plot(z2,power,'-')
plt.xlabel(r'$\bar{z}_2$')
plt.ylabel("power (w)")


# stp.savefig(fname+".png")


#%% current  
# plt.figure()

x1 = np.min(z2)
x2 = np.max(z2)
dx = (x2-x1)/hid["beamCurrentSI"].shape[0]
xx = np.arange(x1,x2,dx)

ax2 = ax1.twinx()
ax2.plot(xx,hid["beamCurrentSI"][:],'--r')
plt.xlabel(r'$\bar{z}_2$')
plt.ylabel("current (A)")

stp.savefig(fname+"_current_power.png")


#%% plot current & power evolution

jjl = np.arange(0,121,10)

for jj in jjl:
    # fname = "fig7_integrated_{jj}.h5"
    fname = f"lcls_integrated_{jj}.h5"
    
    
    hid = h5py.File(fname)
    
    h5f = tables.open_file(fname, mode='r')
    zD = h5f.root.runInfo._v_attrs.zTotal
    
    #% power profile in SI unit
    # nz2 =  h5f.root.runInfo._v_attrs.nZ2
    # dz2 =  h5f.root.runInfo._v_attrs.sLengthOfElmZ2
    # lc  =  h5f.root.runInfo._v_attrs.Lc
    
    spos = hid["runInfo"].attrs["zTotal"]
    
    nz2 = hid["runInfo"].attrs["nZ2"]
    dz2 = hid["runInfo"].attrs["sLengthOfElmZ2"]
    lc  = hid["runInfo"].attrs["Lc"]
    
    c0  =  2.99792458e8
    
    lenz2 = (nz2-1) *dz2
    z2 = (np.arange(0,nz2)) *dz2
    s = z2 *lc
    t = s / c0
    
    power = hid["powerSI"][:]
    
    # plt.figure()
    # plt.title(f"s={spos:.2f} (m)")
    # plt.plot(t,power,'-')
    # plt.xlabel('t (s)')
    # plt.ylabel("power (w)")
    
    # plt.figure()
    
    fig, ax1 = plt.subplots()
    
    plt.title(f"s={spos:.2f} (m)")
    ax1.plot(z2,power,'-')
    plt.xlabel(r'$\bar{z}_2$')
    plt.ylabel("power (w)")
    
    #% current  
    x1 = np.min(z2)
    x2 = np.max(z2)
    dx = (x2-x1)/hid["beamCurrentSI"].shape[0]
    xx = np.arange(x1,x2,dx)
    
    ax2 = ax1.twinx()
    ax2.plot(xx,hid["beamCurrentSI"][:],'--r')
    plt.xlabel(r'$\bar{z}_2$')
    plt.ylabel("current (A)")
    
    stp.savefig(fname+"_current_power.png")


#%% gain curve
import h5py 

hid = h5py.File("lcls_integrated_all.vsh5")
# hid = h5py.File("lcls_integrated_all.vsh5")

plt.figure()
plt.plot(hid["PeakPower"][:])
# plt.plot(hid["power_SI"][:])
# plt.plot(hid["power_SI_Norm"][:])
# plt.plot(hid["Energy"][:])
plt.yscale("log")

