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
    
    field_maps = os.path.join(rootdir, 'sync', 'field-maps')
    gun_profile = 'gun51cavity.txt'
    
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
    print('There is no interpolation function for gun and booster!')


astra_to_sco = 1.590444459638447
TDS_1MV_ratio = 0.07345637 # TO use this ratio, the input Efield will be the TDS voltage in the unit of MV
ast_fmt = '%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%4d%4d'

# # PITZ coordinates
# Lquad = 0.0675
# quads = {'HIGH1.Q1':4.79, 'HIGH1.Q2':5.005, 'HIGH1.Q3':5.6025, 'HIGH1.Q4':5.8525, 'HIGH1.Q5':6.6475,
#         'HIGH1.Q6':6.8925, 'HIGH1.Q7':8.18, 'HIGH1.Q8':8.655, 'HIGH1.Q9':10.208, 'HIGH1.Q10':10.388,
#         'PST.QM1':12.088, 'PST.QM2':12.468, 'PST.QM3':12.848, 'PST.QT1':13.228, 'PST.QT2':13.608,
#         'PST.QT3':13.988, 'PST.QT4':14.368, 'PST.QT5':14.748, 'PST.QT6':15.128, 'HIGH2.Q1':16.635,
#         'HIGH2.Q2':16.735, 'HIGH2.Q3':19.587, 'HIGH2.Q4':21.600, 'HIGH2.Q5':23.013, 'HIGH3.Q1':26.350,
#         'HIGH3.Q2':26.750, 'HIGH3.Q3':27.150} # 2019.10.05

# quads.update({'HIGH2.Q3':19.587, 'HIGH2.Q4':21.600, 'HIGH2.Q5':23.013,
#               'HIGH3.Q1':26.350, 'HIGH3.Q2':26.750, 'HIGH3.Q3':27.150})

# quads.update({'BIO.Q1':28.10, 'BIO.Q2':28.25, 'BIO.Q3':28.40, 
#               'V.Q1':0.23, 'V.Q2':1.5175, 'V.Q3':2.5575, 
#               'TEMP.Q1':5.25, 'TEMP.Q2':5.50, 'TEMP.Q3':5.75}) # 2020

# quads.update({'HIGH2.Q3':19.587, 'HIGH2.Q4':22.1935, 'HIGH2.Q5':23.0785, 
#               'HIGH3.Q1':27.0185, 'HIGH3.Q2':27.4185, 'HIGH3.Q3':27.8185,
#               'BACK3.Q1':0.2685, 'BACK3.Q2':0.2685+0.4, 'BACK3.Q3':0.2685+0.4*2}) # 2021.04.21

# quads.update({'HIGH2.Q3':19.250, 'HIGH2.Q4':21.460, 'HIGH2.Q5':23.220, 
#               'HIGH3.Q1':27.108, 'HIGH3.Q2':27.338, 'HIGH3.Q3':27.778}) # 2021.08.17

# # 2022.01.04, for backward simulation only
# quads['BACK.Q1'] = 29.887-1.8-quads['HIGH3.Q3']
# quads['BACK.Q2'] = quads['HIGH3.Q3']-quads['HIGH3.Q2']+quads['BACK.Q1']
# quads['BACK.Q3'] = quads['HIGH3.Q2']-quads['HIGH3.Q1']+quads['BACK.Q2']

# # 2022.02.02, for backward simulation only
# zmatch = 25.293
# quads['BACK.QM1'] = zmatch - quads['HIGH2.Q5']
# quads['BACK.QM2'] = zmatch - quads['HIGH2.Q4']
# quads['BACK.QM3'] = zmatch - quads['HIGH2.Q3']
# quads['BACK.QM4'] = zmatch - quads['HIGH2.Q2']

# # 2021.10.27, for simulation only
# quads['BIO.QM0'] = 0.1
# quads['BIO.QM1'] = 0.2
# quads['BIO.QM2'] = 0.5
# quads['BIO.QM3'] = 0.8

# quads['BIO.QM4'] = 2.05
# quads['BIO.QM5'] = 2.35
# quads['BIO.QM6'] = 2.65

# # 2023.04.27, for simulation only
# quads['BIO.Q1'] = 0.340
# quads['BIO.Q2'] = 0.655
# quads['BIO.Q3'] = 1.105
# quads['BIO.Q4'] = 1.420


# scrns = {'LOW.SCR1':0.8030, 'LOW.SCR2':1.3790, 'LOW.SCR3':1.7080, 'HIGH1.SCR1':5.2770, 'HIGH1.SCR2':6.2500,\
#         'HIGH1.SCR3':7.1250, 'HIGH1.SCR4':8.4100, 'HIGH1.SCR5':8.9200, 'PST.SCR1':12.2780, 'PST.SCR2':13.0380,\
#         'PST.SCR3':13.7980, 'PST.SCR4':14.5580, 'PST.SCR5':15.3180, 'HIGH2.SCR1':16.3030, 'HIGH2.SCR2':18.2620,\
#         'HIGH2.SCR3':22.86, 'HIGH3.SCR1':26.950, 'BACK2.SCR1':0.46, 'BACK3.SCR1':5.227}
    
# scrns = {'LOW.SCR1':0.8030, 'LOW.SCR2':1.3790, 'LOW.SCR3':1.7080, 'HIGH1.SCR1':5.2770, 'HIGH1.SCR2':6.2500,\
#         'HIGH1.SCR3':7.1250, 'HIGH1.SCR4':8.4100, 'HIGH1.SCR5':8.9200, 'PST.SCR1':12.2780, 'PST.SCR2':13.0380,\
#         'PST.SCR3':13.7980, 'PST.SCR4':14.5580, 'PST.SCR5':15.3180, 'HIGH2.SCR1':16.3030, 'HIGH2.SCR2':18.2620,\
#         'HIGH2.SCR3':23.450, 'HIGH3.SCR1':27.627, 'HIGH3.UND':28.087, 'BACK2.SCR1':0.46, 'BACK3.SCR1':5.227}

# # 2022.07.14
# scrns['HIGH3.SCR1'] = 27.558
# scrns['HIGH3.SCR2'] = 32.040
# scrns['HIGH3.SCR3'] = 33.040

# # 2021.10.27
# scrns['BIO.SCR1'] = 1.425
# scrns['BIO.WINDOW'] = 4.275

# # 2021.11.07
# scrns['HIGH3.D1_ENT'] = 25.793
# scrns['HIGH3.D2_EXIT'] = 27.4269

# # Imaging quads after exit window
# Ldrift = 0.12
# quads['BIO.QI1'] = scrns['BIO.WINDOW']+Ldrift*1+Lquad/2.0
# quads['BIO.QI2'] = scrns['BIO.WINDOW']+Ldrift*2+Lquad*1.0+Lquad/2.0
# quads['BIO.QI3'] = scrns['BIO.WINDOW']+Ldrift*4+Lquad*2.0+Lquad/2.0
# quads['BIO.QI4'] = scrns['BIO.WINDOW']+Ldrift*5+Lquad*3.0+Lquad/2.0

# scrns['BIO.WATER'] = scrns['BIO.WINDOW']+Ldrift*6+Lquad*4

# scrns['BACK.SCR1'] = quads['BACK.Q3']+(quads['HIGH3.Q1']-scrns['HIGH3.SCR1'])
# scrns['BACK.SCR2'] = quads['BACK.Q3']+(quads['HIGH3.Q1']-scrns['HIGH2.SCR3'])

# # quads from the exit of High3.D1, the first dipole of the dogleg
# l1, l2, l3 = 0.45, 0.379286, 0.15
# quads['DOGLEG.Q1'] = l1+Lquad/2
# quads['DOGLEG.Q2'] = quads['DOGLEG.Q1']+l2+Lquad
# quads['DOGLEG.Q3'] = quads['DOGLEG.Q2']+l3*2+Lquad
# quads['DOGLEG.Q4'] = quads['DOGLEG.Q3']+l2+Lquad
# scrns['DOGLEG.LEND'] = 2*(l1+Lquad+l2+Lquad+l3)

# # quads backwards from around window
# quads['BACK.QW1'] = 0.05
# quads['BACK.QW2'] = quads['BACK.QW1']+0.1
# quads['BACK.QW3'] = quads['BACK.QW2']+0.1


# # 2024.04.15
# scrns['DISP5.SCR1'] = 30.109 # length of path from cathode

# steerers = {'LOW.ST1':0.4920, 'LOW.ST2':0.9630, 'LOW.ST3':1.2700, 
#             'LOW.ST4':1.4630, 'LOW.ST5':1.9860, 'LOW.ST5-A':1.9860, 'LOW-HORI.ST5':1.9860, 
#             'LOW.ST5-B':1.9860, 'LOW-VERT.ST5':1.9860, 
#             'HIGH1.ST1':4.8950, 'HIGH1.STA1':5.4270, 'HIGH1.ST2':7.0400, 'HIGH1.STA2':7.2975, 
#             'HIGH1.ST3':8.1280, 'HIGH1.ST4':9.8230, 
#             'PST.ST1':11.600, 'PST.ST2':12.390, 'PST.ST3':13.150, 
#             'PST.ST4':13.910, 'PST.ST5':14.670, 'PST.ST6':14.981, 
#             'HIGH2.ST1':16.453, 'HIGH2.ST2':16.832, 'HIGH2.ST3':19.120, 'HIGH2.ST5':23.042, 
#             'HIGH3.ST1':25.685, 'HIGH3.ST2':27, 'HIGH3.ST3':27.240}
# # 2024.04.15
# steerers['DISP5.D1'] = 26.257 # np.pi/3/60*18*1000

# rotational = ['LOW.ST1', 'LOW.ST3', 'LOW.ST4', 'LOW.ST5', 'HIGH2.ST1', 'HIGH2.ST3', 'HIGH2.ST5', 'HIGH3.ST1', 'HIGH3.ST2', 'HIGH3.ST3']
# for name in rotational:
#     steerers[name+'.IX'] = steerers[name]
#     steerers[name+'.IY'] = steerers[name]
    
# bpms = {'LOBPM1':0.664, 'LOBPM2':1.319, 'BOOBPM1':2.455,'BOOBPM2':4.605,
#       'RFDBPM1':10.599,'RFDBPM2':11.3718,'PSTBPM1':12.2013,'PSTBPM2':12.9613,
#       'PSTBPM3':13.7213,'PSTBPM4':14.4813,'PSTBPM5':15.2413,'HI2BPM1':22.889,
#       'HI3BPM1':26.834,'HI3BPM2':28.029,'HI3BPM3':31.745}


# # 2024.04.15
# # 2024.04.16 commented them to have better trajectory picture during the transport. 
# # TODO: We need to implement dome kind of jumper from one setup to another
# # Just managed. Now one has to change the labels accordingly at trajectory plot
# bpms['DISP5.BPM1'] = 27.407 # Length of path from cathode
# bpms['DISP5.BPM2'] = 29.831 # Length of path from cathode


# bpms.update({'LOW.BPM1':0.664, 'LOW.BPM2':1.319, 'BOOSTER.BPM1':2.455,'BOOSTER.BPM2':4.605,
#       'RFD.BPM1':10.599,'RFD.BPM2':11.3718,'PST.BPM1':12.2013,'PST.BPM2':12.9613,
#       'PST.BPM3':13.7213,'PST.BPM4':14.4813,'PST.BPM5':15.2413,'HIGH2.BPM1':22.889,
#       'HIGH3.BPM1':26.834,'HIGH3.BPM2':28.029,'HIGH3.BPM3':31.745})

# try: # Redefine as case insensitive dictionary
#     from requests.structures import CaseInsensitiveDict
#     pitz = CaseInsensitiveDict(**scrns, **quads, **steerers, **bpms)
#     scrns = CaseInsensitiveDict(**scrns)
#     quads = CaseInsensitiveDict(**quads)
#     steerers = CaseInsensitiveDict(**steerers)
#     bpms = CaseInsensitiveDict(**bpms)
    
# except Exception as err:
#     print(err)
#     pitz = {}
#     pitz.update(**scrns, **quads, **steerers, **bpms)

