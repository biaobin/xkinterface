# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 00:10:39 2021

@author: lixiangk
"""

# from IPython.display import Image
import numpy as np
import scipy as sp
from scipy.interpolate import interp1d, interp2d
from scipy.optimize import curve_fit

from scipy.special import jn, jn_zeros, j0, j1
from scipy.integrate import quad, ode

from scipy.constants import codata

import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.ticker import AutoMinorLocator

import os, re

from timeit import default_timer
import time

from interface import *

def create_folder(y):
    direc = ''
    for i, xi in enumerate(y):
        direc += str.format('Q%d-%.2fT_m-' %  (i+1, xi))
    direc = direc[:-1]
    return direc

def polyfunc(x, *popt):
    x = np.atleast_1d(x)
    r = np.zeros(len(x))
    for i, v in enumerate(popt):
        r += v*x**i
    return r

def polyfit(x, y, *args, **kwargs):
    
    def func(x, a, b, c, d):
        return a+b*x+c*x*x+d*x*x*x
        
    popt, pcov = curve_fit(func, x, y, *args, **kwargs)
    return popt, pcov


#%% THz matching
#%%%   Backward tracking from undultor entrance to High2.Scr3
workdir = r'Backward-3'
os.chdir(workdir)

res = np.loadtxt('Backward-matching-scan.txt')

ix = 3
iy = ix+1

reset_margin()
fig, ax = plt.subplots(); ax.grid()

temp = []
for i, v in enumerate(res[:]):
    #print(v[0:3][::-1])
    direc = create_folder (v[0:3])
    print(direc)

    quads['BACK.Q3'] = 28.09-quads['HIGH3.Q1']
    data = np.loadtxt(direc+os.sep+'BeamDynamics.dat')
    select = data[:,0]>quads['BACK.Q3']+0.0675
    
    data = data[select]
    s1, s2 = scrns['BACK.SCR2'], scrns['BACK.SCR1']
    
    #select = data[:,0]>s2-0.5
    #data = data[select]
    
    popt, pcov = polyfit(data[:,0], np.sqrt(data[:,ix]*data[:,iy]))
    x1, x2 = polyfunc([s1, s2], *popt)
    popt, pcov = polyfit(data[:,0], data[:,iy])
    y1, y2 = polyfunc([s1, s2], *popt)#polyfunc(5.227), polyfunc(0.46)
    
    zz = np.linspace(0, 7)
    ax.plot(data[:,0], data[:,ix], '*')
    ax.plot(zz, polyfunc(zz, *popt), '--')
    
    temp.append([v[0], v[1], v[2], x1, x2])
ax.set_xlabel('RMS size (mm)')
ax.set_ylabel(r'$z$ (m)')
fig.savefig('Backward-tracking.png')

temp = np.array(temp); temp1 = temp
np.savetxt('RMS-size-with-focusing.dat', temp, fmt = '%12.6f')

fig, ax = plt.subplots(figsize = (4, 4))
ax.plot(temp[:,3], temp[:,4], '-o')
ax.grid()
ax.set_xlabel(r'RMS size at High2.Scr3 (mm)')
ax.set_ylabel(r'RMS size at High3.Scr1 (mm)')

fig.savefig('RMS-size-with-focusing.png')

os.chdir('..')

#%%% Forward tracking from EMSY3 to High3.Scr1
workdir = r'Forward-T3'
os.chdir(workdir)

res = np.loadtxt('Forward-focusing-scan.txt')
z0 = 18.262
#z0 = 16.3

fig1, ax1 = plt.subplots(); ax.grid()

temp = []
for i, v in enumerate(res[:]):
    #print(v[0:3][::-1])
    direc = create_folder(v)
    print(direc)
    
    data = np.loadtxt(direc+os.sep+'BeamDynamics.dat')
    data[:,0] += z0
    
    fx = interp1d(data[:,0], data[:,ix])
    fy = interp1d(data[:,0], data[:,iy])
    
    ax1.plot(data[:,0], data[:,3], '--')
    ax1.plot(data[:,0], data[:,4], '--')
    s1, s2 = scrns['HIGH2.SCR3'], scrns['HIGH3.SCR1']
    #s1, s2 = 16.3-1, 16.3+1
    temp.append([v[0], v[1], v[2], np.sqrt(fx(s1)*fy(s1)), np.sqrt(fx(s2)*fy(s2))])
ax1.set_xlabel('RMS size (mm)')
ax1.set_ylabel(r'$z$ (m)')
fig1.savefig('Forward-tracking.png')


temp = np.array(temp); temp2 = temp
np.savetxt('RMS-size-with-focusing.dat', temp, fmt = '%12.6f')

#fig, ax = plt.subplots(figsize = (5, 4))
ax.plot(temp[:,3], temp[:,4], '-o')
#ax.grid()

#ax.set_xlim(0.5, 2.5)
#ax.set_ylim(0.9, 1.5)
ax.legend(['Backward tracking', 'Forward tracking'], loc = 'upper left')

if ix == 3:
    ax.set_xlabel(r'RMS size at Scr1 (mm)')
    ax.set_ylabel(r'RMS size at Scr2 (mm)')
    fig.savefig('RMS-size-with-focusing2.png')
elif ix == 6:
    ax.set_xlabel(r'$\beta$ function at High2.Scr3 (m)')
    ax.set_ylabel(r'$\beta$ function at High3.Scr1 (m)')
    fig.savefig('beta-function-with-focusing2.png')

os.chdir('..')

#%%% Get the intersection of the two curves
from intersect import intersection

x, y = intersection(temp1[:,3], temp1[:,4], temp2[:,3], temp2[:,4])

T2 = np.zeros(3)
for i in [0, 1, 2]:
    fQ = interp1d(temp2[:,3], temp2[:,i], kind = 'quadratic')
    T2[i] = fQ(x)[0]
print('Forward triplet: ', T2)

T1 = np.zeros(3)
for i in [0, 1, 2]:
    fQ = interp1d(temp1[:,3], temp1[:,i], kind = 'quadratic')
    T1[i] = fQ(x)[0]
print('Backward triplet in forward order: ', T1[::-1])
