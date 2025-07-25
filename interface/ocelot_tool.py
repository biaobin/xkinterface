# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 09:51:10 2022

@author: lixiangk
"""
# the output of plotting commands is displayed inline within frontends, 
# directly below the code cell that produced it
import time 

# Ocelot packages
# this python library provides generic shallow (copy) and deep copy (deepcopy) operations 
from copy import deepcopy
# import from Ocelot main modules and functions
from ocelot import *
# import from Ocelot graphical modules
# from ocelot.gui.accelerator import *
# load beam distribution
# this function convert CSRtrack beam distribution to Ocelot format
# - ParticleArray. ParticleArray is designed for tracking.
# in order to work with converters we have to import 
# specific module from ocelot.adaptors
from ocelot.adaptors.csrtrack2ocelot import *
from ocelot.adaptors.astra2ocelot import *

def G2K(G, P0, qn = 1.):
    '''
    Parameters
      G: gradient, T/m
      P0: momentum, MeV/c
      qn = q/qe: number of charge
    Returns
      K: focal strength, m^-2
    '''
    ccc = 299792458.0
    return G*qn/(P0*1e6/ccc)

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_twiss(tws_track, z0 = 0, nrows = 1, label = '', fig_ext = None, fout = None, plot = True):
    
    bg = (tws_track[0].E*1e3)/0.51099895000; print(bg)
        
    res = []
    for i, tw in enumerate(tws_track):
        res.append([tw.s, tw.emit_x*bg*1e6, tw.emit_y*bg*1e6, 
                    tw.xx**0.5*1e3, tw.yy**0.5*1e3, tw.tautau**0.5*1e3,
                    tw.beta_x, tw.beta_y, tw.alpha_x, tw.alpha_y])
    res = np.array(res)
    res[:,0] += z0
    
    if fout != None:
        fout = 'track-%s.dat' % label
        np.savetxt(fout, res, fmt = '%14.6E')
    
    if plot:
        fig, ax = plt.subplots(nrows = nrows, ncols = 4//nrows, figsize = (6*4/nrows, 6*nrows))
        ax = ax.flatten()
        
        ax[0].plot(res[:,0], res[:,1], 'r-', label = r'$x$')
        ax[0].plot(res[:,0], res[:,2], 'b-', label = r'$y$')
        
        
        
        ax[1].plot(res[:,0], res[:,3], 'r-', label = r'$x$')
        ax[1].plot(res[:,0], res[:,4], 'b-', label = r'$y$')
        ax[1].plot(res[:,0], res[:,5], 'g-', label = r'$z$')
        
        #ax[1].plot(res[:,0], np.sqrt(res[:,3]*res[:,4]), 'k-o', label = r'$xy$')
        
        ax[2].plot(res[:,0], res[:,6], 'r-', label = r'$x$')
        ax[2].plot(res[:,0], res[:,7], 'b-', label = r'$y$')
        
        ax[3].plot(res[:,0], res[:,8], 'r-', label = r'$x$')
        ax[3].plot(res[:,0], res[:,9], 'b-', label = r'$y$')
        
        ylabels = ['Norm. emit. (um)', 'RMS size (mm)', 'beta function (m)', 'alpha function']
        titles = ['Emittance', 'RMS size', 'beta', 'alpha']
        for i in np.arange(len(ax)):
            ax[i].set_title(titles[i] + ' vs s')
            ax[i].set_xlabel(r'$s$ (m)')
            ax[i].set_ylabel(ylabels[i])
            ax[i].legend()
            ax[i].grid(True)
        plt.tight_layout()
        # plt.grid(True)
        if fig_ext != None:
            fig.savefig(('%s' % label)+fig_ext)
    return res

def plot_sig(tws_track, z0 = 0, nrows = 1, label = '', fig_ext = None, fout = None, plot = True):
    
    bg = (tws_track[0].E*1e3)/0.51099895000; print(bg)
        
    res = []
    for i, tw in enumerate(tws_track):
        res.append([tw.s, tw.emit_x*bg*1e6, tw.emit_y*bg*1e6, 
                    tw.xx**0.5*1e3, tw.yy**0.5*1e3, tw.tautau**0.5*1e3,
                    tw.beta_x, tw.beta_y, tw.alpha_x, tw.alpha_y])
    res = np.array(res)
    res[:,0] += z0
    
    if fout != None:
        fout = 'track-%s.dat' % label
        np.savetxt(fout, res, fmt = '%14.6E')
    
    if plot:
                
        # plt.figure(figsize=(8,6))
        plt.plot(res[:,0], res[:,3], 'r-', label = r'ocelot, sigx')
        plt.plot(res[:,0], res[:,4], 'b-', label = r'ocelot, sigy')
        # plt.plot(res[:,0], res[:,5], 'g-', label = r'ocelot, sigz')
        
        plt.xlabel("s (m)")
        plt.ylabel("RMS size (mm)")
        plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.grid(True)
        if fig_ext !=None:
            plt.savefig(('%s' % label)+fig_ext)   #fig_ext -> ".png"
    return res

def plot_result(res, z0 = 0, nrows = 1, label = '', fig_ext = '.png', fout = None):
    
    res = np.array(res)
    #res[:,0] += z0
    
    if fout != None:
        fout = 'track-%s.dat' % label
        np.savetxt(fout, res, fmt = '%14.6E')
    
    fig, ax = plt.subplots(nrows = nrows, ncols = 4//nrows, figsize = (3*4/nrows, 3*nrows))
    ax = ax.flatten()
    
    ax[0].plot(res[:,0], res[:,1], 'r-*', label = r'$x$')
    ax[0].plot(res[:,0], res[:,2], 'b-<', label = r'$y$')
    
    ax[1].plot(res[:,0], res[:,3], 'r-*', label = r'$x$')
    ax[1].plot(res[:,0], res[:,4], 'b-<', label = r'$y$')
    ax[1].plot(res[:,0], res[:,5], 'g-<', label = r'$z$')
    
    ax[2].plot(res[:,0], res[:,6], 'r-*', label = r'$x$')
    ax[2].plot(res[:,0], res[:,7], 'b-<', label = r'$y$')
    
    ax[3].plot(res[:,0], res[:,8], 'r-*', label = r'$x$')
    ax[3].plot(res[:,0], res[:,9], 'b-<', label = r'$y$')
    
    ylabels = ['Norm. emit. (um)', 'RMS size (mm)', 'beta function (m)', 'alpha function']
    titles = ['Emittance', 'RMS size', 'beta', 'alpha']
    for i in np.arange(len(ax)):
        ax[i].set_title(titles[i] + ' vs s')
        ax[i].set_xlabel(r'$s$ (m)')
        ax[i].set_ylabel(ylabels[i])
        ax[i].legend()
        ax[i].grid()
    plt.tight_layout()
    
    fig.savefig(('%s' % label)+fig_ext)
    return res

def plot_particleArray(p_array, z0 = 0, xmax = 5, ncols = 5, nrows = 5, fig_ext = '.png', **kwargs):
    
    nn = len(p_array)
    step = nn//(ncols*nrows); print(nn, step)
    if step < 1:
        step = 1
    #if step == 1:
    #    step = 2
    
    fig, ax = plt.subplots(ncols = ncols, nrows = nrows, figsize = (3*nrows, 3*ncols), **kwargs)
    ax = ax.flatten()
    
    #step = 8
    # xmax = 8
    # if 'xmax' in kwargs.keys():
    #     xmax = kwargs['xmax']
    # else:
    #     xmax = 5
    ymax = xmax
    for i in np.arange(0, len(p_array), step):
        k = i//step
        if k >= ncols*nrows:
            print('Subplot indexing out of range!')
            k = -1
            
        if i+step >= len(p_array):
            j = -1
        else:
            j = i
            
        ax[k].hist2d(p_array[j].x()*1000, p_array[j].y()*1000, bins = 200,
                           cmin = 1, range= [[-xmax, xmax], [-ymax, ymax]])
        ax[k].set_title('s = %.2f m' % (p_array[j].s+z0))
        ax[k].set_xlim(-xmax, xmax)
        ax[k].set_ylim(-ymax, ymax)
        
        ax[k].set_aspect('equal')
        if j == -1:
            #xmax = 5
            ax[k].set_xlim(-xmax, xmax)
            ax[k].set_ylim(-xmax, xmax)   
            
    plt.tight_layout()
    fig.savefig('xy2d-along-beamline'+fig_ext)
    return

def plot_particleArray_zpz(p_array, z0 = 0, xmax = 5, ncols = 5, nrows = 5, fig_ext = '.png', **kwargs):
    
    nn = len(p_array)
    step = nn//(ncols*nrows); print(nn, step)
    #if step == 1:
    #    step = 2
    
    fig, ax = plt.subplots(ncols = ncols, nrows = nrows, figsize = (3*nrows, 3*ncols), **kwargs)
    ax = ax.flatten()
    
    #step = 8
    # xmax = 8
    # if 'xmax' in kwargs.keys():
    #     xmax = kwargs['xmax']
    # else:
    #     xmax = 5
    ymax = xmax
    for i in np.arange(0, len(p_array), step):
        k = i//step
        if k >= ncols*nrows:
            print('Subplot indexing out of range!')
            k = -1
            
        if i+step >= len(p_array):
            j = -1
        else:
            j = i
            
        ax[k].hist2d(p_array[j].tau()*1000, p_array[j].p()*1000, bins = 200,
                           cmin = 1, range= [[-xmax, xmax], [-ymax, ymax]])
        ax[k].set_title('s = %.2f m' % (p_array[j].s+z0))
        #ax[k].set_xlim(-xmax, xmax)
        #ax[k].set_ylim(-ymax, ymax)
        
        #ax[k].set_aspect('equal')
        #if j == -1:
        #    #xmax = 5
        #    ax[k].set_xlim(-xmax, xmax)
        #    ax[k].set_ylim(-xmax, xmax)   
            
    plt.tight_layout()
    fig.savefig('zpz-along-beamline'+fig_ext)
    return
