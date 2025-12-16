import numpy as np
import matplotlib.pyplot as plt

def getchirp(p_array):
    # return chirp h [/m]
    x = -p_array.tau()
    y = p_array.p()

    coeff = np.polyfit(x, y, deg=5)
    
    return coeff[4]

def removechirp(x,y,removeorder=2):
    # remove 1,2,...removeorder energy modulation
    coeff = np.polyfit(x, y, deg=5)   # a,b,c
    
    #for order<= removeto, remove it
    coeff2 = coeff[5-removeorder:]
    poly = np.poly1d(coeff2)
    
    y_detrend = y - poly(x)
    
    return x,y_detrend

def longi_show(p_array):
    # with chirp h labeled
    #-------------------------
    plt.figure(figsize=(15,6))
    
    plt.subplot(1,2,1)
    # plt.title(f"theta={angle/np.pi*180:.1f} deg")
    x = -p_array.tau()
    y = p_array.p()

    coeff = np.polyfit(x, y, deg=5)
    
    plt.scatter(x/2.998e8*1e12, y, s=0.01, color='r',label=f"h={coeff[4]:.2f}/m")
    # impzplt.plot_heatpha(x, y_detrend, bins=100)
    
    plt.xlabel('z (ps)')
    plt.ylabel('$\Delta E/E$')
    plt.legend(markerscale=40)

    plt.subplot(1,2,2)
    plt.plot(-p_array.I()[:,0]/2.998e8*1e12, p_array.I()[:,1], 'r-')
    plt.xlabel('z (ps)')
    plt.ylabel('current (A)')

def compare_longi_array12(p_array_i, p_array,label1="before BC",label2="after BC"):
    plt.figure(figsize=(15,6))
    
    plt.subplot(1,2,1)
    # plt.title(f"theta={angle/np.pi*180:.1f} deg")
    # plt.plot(-p_array_i.tau()/2.998e8 *1e12, p_array_i.p(), 'go')
    # plt.plot(-p_array.tau()/2.998e8 *1e12, p_array.p(), 'r.')
    
    plt.scatter(-p_array_i.tau()/2.998e8 *1e12, p_array_i.p(), color='g',s=0.01,label=label1)
    plt.scatter(-p_array.tau()/2.998e8 *1e12, p_array.p(),     color='r',s=0.01,label=label2)
    
    # impzplt.plot_heatpha(-p_array_i.tau()/2.998e8 *1e12, p_array_i.p(), bins=100)
    # impzplt.plot_heatpha(-p_array.tau()/2.998e8 *1e12, p_array.p(), bins=100)
    
    plt.xlabel('z (ps)')
    plt.ylabel('$\Delta E/E$')
    plt.legend(markerscale=40)
    
    plt.subplot(1,2,2)
    plt.plot(-p_array_i.I()[:,0]/2.998e8*1e12, p_array_i.I()[:,1], 'g--',label=label1)
    plt.plot(-p_array.I()[:,0]/2.998e8*1e12, p_array.I()[:,1], 'r-',label=label2)
    plt.xlabel('z (ps)')
    plt.ylabel('current (A)')
    # plt.gca().invert_xaxis()
    
    plt.legend(markerscale=40)
    
    plt.show() 
    plt.grid(True)
    plt.show()

def longi_removechirp(p_array, removeorder=1):
    # ------------------------------------------------------------
    plt.figure(figsize=(15,6))
    
    plt.subplot(1,2,1)
    # plt.title(f"theta={angle/np.pi*180:.1f} deg")
    x = -p_array.tau()/2.998e8 *1e12
    y = p_array.p()
    x, y_detrend = removechirp(x, y, removeorder=removeorder)
    
    plt.scatter(x, y, s=0.01, color='g',label='initial')
    plt.scatter(x, y_detrend, s=0.01, color='r',label='chirp removed')
    # impzplt.plot_heatpha(x, y_detrend, bins=100)
    
    plt.xlabel('z (ps)')
    plt.ylabel('$\Delta E/E$')
    plt.legend(markerscale=40)

    plt.subplot(1,2,2)
    plt.plot(-p_array.I()[:,0]/2.998e8*1e12, p_array.I()[:,1], 'r-')
    plt.xlabel('z (ps)')
    plt.ylabel('current (A)')
    
    plt.grid(True)
    plt.show()

from ocelot import * 
def get_BC_r56(cell, energy=17e-3):
    
    lat_chic = MagneticLattice(cell)
    # in that case energy is not important we do not have 
    # energy dependant elements here
    R = lattice_transfer_map(lat_chic, energy=energy)
    print("\nR56 = ", R[4,5]*1000, "mm")
