import sys

from timeit import default_timer
import time

from scipy.interpolate import interp1d, interp2d
from scipy.optimize import curve_fit
import numpy as np
import os

# Load coordinates
from .PITZcoordinates import *

I2B = lambda I: -(0.0000372+0.000588*I)
B2I = lambda B: (-B-0.0000372)/0.000588

try:
    ''' 
    Load gun and booster 2D (grad and phase) data to create interpolation 
    functions for gun and booster
    - get_MaxE_gun(phi_gun = 0, EG = 6.3) 
      returns gun gradient for given gun phase and beam momentum
    - get_MaxE_booster(MaxE_gun = 60, phi_gun = 0, phi_booster = 0, P = 17.05)
      returns booster gradient for given gun grad&phase, booster phase and beam momentum
    '''
    
    cwd = os.path.dirname(os.path.abspath(__file__))

    field_maps = os.path.join(cwd, 'field-maps')
    gun_profile = 'gun51cavity.txt'
    boo_profile = 'CDS14_15mm.txt'
    sol_profile = 'gunsolenoidsPITZ.txt'
    
    ### 2D gun gradient and phase scan and interpolation, updated 27.02.2023
    data_gun = np.loadtxt(field_maps+os.sep+'gun51_scan2d.dat')
    shape = (51, 91)
    E_gun   = np.reshape(data_gun[:,0], shape)
    phi_gun = np.reshape(data_gun[:,1], shape)
    EG_gun  = np.reshape(data_gun[:,2], shape)

    # fEG2d_gun = interp2d(E_gun, phi_gun, EG_gun, bounds_error = False)
    from scipy.interpolate import RegularGridInterpolator
    fEG2d_gun_t = RegularGridInterpolator((np.unique(E_gun), np.unique(phi_gun)), EG_gun,
                                    bounds_error = False)
    fEG2d_gun = lambda E, phi: fEG2d_gun_t(([E], [phi]))[0]

    def get_MaxE_gun(phi_gun = 0, EG = 6.3):
        
        EG1 = np.array([fEG2d_gun(E0, phi_gun) for E0 in np.unique(E_gun)])
        fMaxE = interp1d(EG1, np.unique(E_gun),
                         bounds_error = False, fill_value = 'extrapolate')
        
        #EG1 = fEG2d_gun(phi_gun, EG)
        #MaxE_gun = fEgun2d(phi_gun, EG)
        MaxE_gun = fMaxE(EG)
        return MaxE_gun.item()
    ###

    ### 2D booster gradient and phase scan and interpolation, updated 01.03.2023
    # Somehow, the three columns in booster 2d scan data is phi, grad, and momentum,
    # which is different from that for the gun
    data_booster = np.loadtxt(field_maps+os.sep+'phi2_scan.dat')
    shape = (71, 121)
    phi_booster = np.reshape(data_booster[:,0], shape)
    E_booster   = np.reshape(data_booster[:,1], shape)
    EG_booster  = np.reshape(data_booster[:,2], shape)

    # fEG2d_booster = interp2d(E_booster, phi_booster, EG_booster)
    fEG2d_booster_t = RegularGridInterpolator((np.unique(E_booster), np.unique(phi_booster)), EG_booster.T,
                                              bounds_error = False)
    fEG2d_booster = lambda E, phi: fEG2d_booster_t(([E], [phi]))[0]

    def get_MaxE_booster(MaxE_gun = 60, phi_gun = 0, phi_booster = 0, Ef = 17.05):
        
        Eb = np.linspace(10, 22, 121*5)
        EG = fEG2d_booster(Eb, phi_booster)
        fEb_EG = interp1d(EG, Eb)
        
        E1 = fEG2d_gun(MaxE_gun, phi_gun)
        
        MaxE_booster = fEb_EG(Ef-E1)
        return MaxE_booster.item()
    ###
except Exception as err:
    print(err)
    print('There is no interpolation function for gun and booster!')
    
try:
    ### 2D gun gradient and phase scan and interpolation for thermal imaging, updated 07.2024
    data_thermal_imaging = np.loadtxt(field_maps+os.sep+'ThermalImagingScan_gun5.1.dat')
    shape1 = (9, 3)
    E_gun1   = np.reshape(data_thermal_imaging[:,0], shape1)
    phi_gun1 = np.reshape(data_thermal_imaging[:,1], shape1)
    Imain_te  = np.reshape(data_thermal_imaging[:,2], shape1)

    fTE2d_gun_t = RegularGridInterpolator((np.unique(E_gun1), np.unique(phi_gun1)),
                                          Imain_te,
                                          bounds_error = False)
    def fTE2d_gun(E, phi): 
        return fTE2d_gun_t(([E], [phi]))[0]
    ###
    
except Exception as err:
    print(err)
    print('There is no interpolation function for thermal emittance!')


astra_to_sco = 1.590444459638447
TDS_1MV_ratio = 0.07345637 # TO use this ratio, the input Efield will be the TDS voltage in the unit of MV
ast_fmt = '%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%4d%4d'

def CreateFolderName(x, flag = 'injector', **kwargs):
    Ipart = 250000
    if len(kwargs)>0:
        if 'Ipart' in kwargs.keys():
            Ipart = kwargs['Ipart']
            
    if flag.upper() in ['INJECTOR', 'PI']:
        
        sigma_x = sigma_y = x[1]/4.
        phi_gun, phi_booster = x[3], x[5]
        Imain = x[6]

        MaxE_gun = x[2]
        phi_gun, phi_booster = x[3], x[5]
        MaxE_gun = get_MaxE_gun(phi_gun, 6.3)
    
        MaxE_booster = get_MaxE_booster(MaxE_gun, phi_gun, phi_booster, 17)
        
        Q_total = x[0]/1e3
        #Ipart = int(Q_total*50e3)

        direc = str.format('Q-%.2fpC-D-%.2fmm-E1-%.2fMV_m-phi1-%.2fdeg-E2-%.2fMV_m-phi2-%.2fdeg-I-%.2fA' %\
                           (Q_total*1e3, x[1], MaxE_gun, phi_gun, MaxE_booster, phi_booster, Imain))
    
    elif flag.upper() in ['QUADS', 'TRANSPORT']:
        direc = ''
        for i in np.arange(len(x[:])):
            direc += str.format('Q%d-%.2fT_m-' %  (i+1, x[i]))
        direc = direc[:-1]
    
    elif flag.upper() in ['TWISS', 'MATCHING']:        
        direc = str.format('n%.0fk-sig_x-%.2fmm-sig_y-%.2fmm-alp_x-%.2f-alp_y-%.2f-nemit_x-%.2fum-nemit_y-%.2fum' % (Ipart/1000., *x))
    
    return direc

def CreateFolder(direc):
    if not os.path.exists(direc):
        os.mkdir(direc)
        print('Create folder: ', direc)
    else:
        print('Folder exits: ', direc)
        
    return