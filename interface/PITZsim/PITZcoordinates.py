# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 17:16:47 2023

@author: lixiangk
"""

import os
import numpy as np

Lquad = 0.0675
quads = {'HIGH1.Q1':4.79, 'HIGH1.Q2':5.005, 'HIGH1.Q3':5.6025, 'HIGH1.Q4':5.8525, 'HIGH1.Q5':6.6475,
        'HIGH1.Q6':6.8925, 'HIGH1.Q7':8.18, 'HIGH1.Q8':8.655, 'HIGH1.Q9':10.208, 'HIGH1.Q10':10.388,
        'PST.QM1':12.088, 'PST.QM2':12.468, 'PST.QM3':12.848, 'PST.QT1':13.228, 'PST.QT2':13.608,
        'PST.QT3':13.988, 'PST.QT4':14.368, 'PST.QT5':14.748, 'PST.QT6':15.128, 'HIGH2.Q1':16.635,
        'HIGH2.Q2':16.735, 'HIGH2.Q3':19.587, 'HIGH2.Q4':21.600, 'HIGH2.Q5':23.013, 'HIGH3.Q1':26.350,
        'HIGH3.Q2':26.750, 'HIGH3.Q3':27.150} # 2019.10.05

quads.update({'HIGH2.Q3':19.250, 'HIGH2.Q4':21.460, 'HIGH2.Q5':23.220, 
              'HIGH3.Q1':27.108, 'HIGH3.Q2':27.338, 'HIGH3.Q3':27.778,
              'HIGH3.Q4':30.108, 'HIGH3.Q5':32.338}) # 2021.08.17
#quads.update({'HIGH3.Q4':30.108, 'HIGH3.Q5':32.338})




scrns = {'LOW.SCR1':0.8030, 'LOW.SCR2':1.3790, 'LOW.SCR3':1.7080, 'HIGH1.SCR1':5.2770, 'HIGH1.SCR2':6.2500,
        'HIGH1.SCR3':7.1250, 'HIGH1.SCR4':8.4100, 'HIGH1.SCR5':8.9200, 'PST.SCR1':12.2780, 'PST.SCR2':13.0380,
        'PST.SCR3':13.7980, 'PST.SCR4':14.5580, 'PST.SCR5':15.3180, 'HIGH2.SCR1':16.3030, 'HIGH2.SCR2':18.2620,
        'HIGH2.SCR3':23.450, 'HIGH3.SCR1':27.627, 'HIGH3.UND':28.087, 'BACK2.SCR1':0.46, 'BACK3.SCR1':5.227}

# 2022.07.14
scrns['HIGH3.SCR1'] = 27.558
scrns['HIGH3.SCR2'] = 32.040
scrns['HIGH3.SCR3'] = 33.040

# 2024.04.15
scrns['DISP5.SCR1'] = 30.109 # length of path from cathode

steerers = {'LOW.ST1':0.4920, 'LOW.ST2':0.9630, 'LOW.ST3':1.2700, 
            'LOW.ST4':1.4630, 'LOW.ST5':1.9860, 'LOW.ST5-A':1.9860, 'LOW-HORI.ST5':1.9860, 
            'LOW.ST5-B':1.9860, 'LOW-VERT.ST5':1.9860, 
            'HIGH1.ST1':4.8950, 'HIGH1.STA1':5.4270, 'HIGH1.ST2':7.0400, 'HIGH1.STA2':7.2975, 
            'HIGH1.ST3':8.1280, 'HIGH1.ST4':9.8230, 
            'PST.ST1':11.600, 'PST.ST2':12.390, 'PST.ST3':13.150, 
            'PST.ST4':13.910, 'PST.ST5':14.670, 'PST.ST6':14.981, 
            'HIGH2.ST1':16.453, 'HIGH2.ST2':16.832, 'HIGH2.ST3':19.120, 'HIGH2.ST5':23.042, 
            'HIGH3.ST1':25.685, 'HIGH3.ST2':27, 'HIGH3.ST3':27.240}
# 2024.04.15
steerers['DISP5.D1'] = 26.257 # np.pi/3/60*18*1000

rotational = ['LOW.ST1', 'LOW.ST3', 'LOW.ST4', 'LOW.ST5', 'HIGH2.ST1', 'HIGH2.ST3', 'HIGH2.ST5', 'HIGH3.ST1', 'HIGH3.ST2', 'HIGH3.ST3']
for name in rotational:
    steerers[name+'.IX'] = steerers[name]
    steerers[name+'.IY'] = steerers[name]
    
bpms = {'LOBPM1':0.664, 'LOBPM2':1.319, 'BOOBPM1':2.455,'BOOBPM2':4.605,
      'RFDBPM1':10.599,'RFDBPM2':11.3718,'PSTBPM1':12.2013,'PSTBPM2':12.9613,
      'PSTBPM3':13.7213,'PSTBPM4':14.4813,'PSTBPM5':15.2413,'HI2BPM1':22.889,
      'HI3BPM1':26.834,'HI3BPM2':28.029,'HI3BPM3':31.745}


# 2024.04.15
# 2024.04.16 commented them to have better trajectory picture during the transport. 
# TODO: We need to implement dome kind of jumper from one setup to another
# Just managed. Now one has to change the labels accordingly at trajectory plot
bpms['DISP5.BPM1'] = 27.407 # Length of path from cathode
bpms['DISP5.BPM2'] = 29.831 # Length of path from cathode

# 2024.11.28
# Quads backward counted from undulator entrance
quads['BACK.Q1'] = scrns['HIGH3.UND']-quads['HIGH3.Q3']
quads['BACK.Q2'] = scrns['HIGH3.UND']-quads['HIGH3.Q2']
quads['BACK.Q3'] = scrns['HIGH3.UND']-quads['HIGH3.Q1']

# Screens backward counted from undulator entrance
scrns['BACK.SCR1'] = scrns['HIGH3.UND']-scrns['HIGH3.SCR1']
scrns['BACK.SCR2'] = scrns['HIGH3.UND']-scrns['HIGH2.SCR3']

bpms.update({'LOW.BPM1':0.664, 'LOW.BPM2':1.319, 'BOOSTER.BPM1':2.455,'BOOSTER.BPM2':4.605,
      'RFD.BPM1':10.599,'RFD.BPM2':11.3718,'PST.BPM1':12.2013,'PST.BPM2':12.9613,
      'PST.BPM3':13.7213,'PST.BPM4':14.4813,'PST.BPM5':15.2413,'HIGH2.BPM1':22.889,
      'HIGH3.BPM1':26.834,'HIGH3.BPM2':28.029,'HIGH3.BPM3':31.745})

zc1 = 18.650
zc2 = 19.630
zc3 = 21.186
zc4 = 22.166
# z-axis coordinate
chicane= {"CHICANE.D1":zc1, "CHICANE.D2":zc2, "CHICANE.D3":zc3, "CHICANE.D4":zc4}

#center of the undulator
undulator= {"UND":29.887}   


try: # Redefine as case insensitive dictionary
    from requests.structures import CaseInsensitiveDict
    pitz = CaseInsensitiveDict(**scrns, **quads, **steerers, **bpms, **chicane, **undulator)
    
    scrns = CaseInsensitiveDict(**scrns)
    quads = CaseInsensitiveDict(**quads)
    steerers = CaseInsensitiveDict(**steerers)
    bpms = CaseInsensitiveDict(**bpms)
    chicane = CaseInsensitiveDict(**chicane)
    undulator = CaseInsensitiveDict(**undulator)
    
except Exception as err:
    print(err)
    pitz = {}
    pitz.update(**scrns, **quads, **steerers, **bpms, **chicane, **undulator)
