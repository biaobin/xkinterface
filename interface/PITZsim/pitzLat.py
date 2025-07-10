from .PITZcoordinates import *
from ocelot import *
import sys
import numpy as np
from collections import defaultdict
from math import *

# to add drift elements automatically between quads, undulator, bend
def add_drift(line: list):
    cell = [ pitzlat[line[0]] ]
    Nelem = len(line)
    for j in range(Nelem-1):
        # NO drift should be added between dipoles before this func is called
        if (pitzlat[line[j]].__class__ in [Bend,SBend]) and (pitzlat[line[j+1]].__class__ in [Bend,SBend]):
            Larc = pitzlat[line[j]].l
            angle = abs(pitzlat[line[j]].angle)
            
            rho = Larc/angle
            Lb = rho*sin(angle)
            e2 = pitzlat[line[j]].e2 
            
            ld = (pitz[line[j+1]] -pitz[line[j]] -Lb)/cos(e2)
        
        elif pitzlat[line[j]].__class__ == Drift:
            print("ERROR, no drift should be added before this func being called.")
            sys.exit()
        else:   
            ld = pitz[line[j+1]] -pitz[line[j]] -0.5*pitzlat[line[j+1]].l -0.5*pitzlat[line[j]].l
        
        cell += [Drift(l=ld), pitzlat[line[j+1]]]   
    return cell

#%% lattice
nest_dict = lambda: defaultdict(nest_dict)
pitzlat = nest_dict()

Lquad, Rquad = 0.0675, 0.0215
astra_to_sco = 1.590444459638447

pitzlat['HIGH1.SCR1'] = Marker()
pitzlat['PST.SCR1']   = Marker()
pitzlat['HIGH2.SCR2'] = Marker()
pitzlat['HIGH2.SCR3'] = Marker()
pitzlat['HIGH3.BPM1'] = Marker()
pitzlat['HIGH3.UND']  = Marker()
pitzlat['HIGH3.SCR2'] = Marker()


pitzlat['HIGH1.Q4'] = Quadrupole(l=Lquad, k1=-0.65/astra_to_sco)
pitzlat['HIGH1.Q6'] = Quadrupole(l=Lquad, k1= 0.91/astra_to_sco)
pitzlat['HIGH1.Q7'] = Quadrupole(l=Lquad, k1=-0.43/astra_to_sco)

pitzlat['PST.QT2'] = Quadrupole(l=Lquad, k1=0.63/astra_to_sco)
pitzlat['PST.QT4'] = Quadrupole(l=Lquad, k1=-1.15/astra_to_sco)
pitzlat['PST.QT6'] = Quadrupole(l=Lquad, k1=0.64/astra_to_sco)

pitzlat['HIGH2.Q2'] = Quadrupole(l=Lquad, k1=0.37/astra_to_sco)
pitzlat['HIGH2.Q5'] = Quadrupole(l=Lquad, k1=-0.3/astra_to_sco)

pitzlat['HIGH3.Q1'] = Quadrupole(l=Lquad, k1=2.30/astra_to_sco)
pitzlat['HIGH3.Q2'] = Quadrupole(l=Lquad, k1=-3.19/astra_to_sco)
pitzlat['HIGH3.Q3'] = Quadrupole(l=Lquad, k1=0.59/astra_to_sco)

# centers of dipoles
Larc = 0.33161
angle = np.pi*19/180
pitzlat["CHICANE.D1"] = Bend(l = Larc, angle=-angle, e1=0.0, e2=-angle, gap=0.08, tilt=pi/2, fint = 0.5, eid='BB.393.B2')
pitzlat["CHICANE.D2"] = Bend(l = Larc, angle= angle, e1=angle,  e2=0.0, gap=0.08, tilt=pi/2, fint = 0.5, eid='BB.402.B2')
pitzlat["CHICANE.D3"] = Bend(l = Larc, angle= angle, e1=0.0,  e2=angle, gap=0.08, tilt=pi/2, fint = 0.5, eid='BB.404.B2')
pitzlat["CHICANE.D4"] = Bend(l = Larc, angle=-angle, e1=-angle, e2=0.0, gap=0.08, tilt=pi/2, fint = 0.5, eid='BB.413.B2')

pitzlat["UND"] = Undulator(3e-2, 113, 3.49)

high1 = ['HIGH1.SCR1', 'HIGH1.Q4', 'HIGH1.Q6', 'HIGH1.Q7', 'PST.SCR1']
pst   = ['PST.SCR1', 'PST.QT2', 'PST.QT4', 'PST.QT6', 'HIGH2.Q2', 'HIGH2.SCR2']
high2 = ['HIGH2.SCR2', 'CHICANE.D1', 'CHICANE.D2', 'CHICANE.D3', 'CHICANE.D4', 'HIGH2.Q5', 'HIGH2.SCR3', 'HIGH3.BPM1']
high3 = ['HIGH3.BPM1', 'HIGH3.Q1', 'HIGH3.Q2', 'HIGH3.Q3', 'HIGH3.UND']
undu  = ['HIGH3.UND', 'UND', 'HIGH3.SCR2']

#pitzlat["high1"] = add_drift(high1)  
#pitzlat["pst"]   = add_drift(pst)
#pitzlat["high2"] = add_drift(high2)
#pitzlat["high3"] = add_drift(high3)

def getLat(line):
    #line should NOT contain drift element
    return add_drift(line)
 
