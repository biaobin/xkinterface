#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:07:20 2025

@author: biaobin
"""
from impzpy.impactz_parser import impactz_parser as impz
from imptpy.impactt_parser import impactt_parser as impt
import os

path='/afs/ifh.de/group/pitz/data/biaobin/sim1/github/pitz_simu/submit_jobs_s2e/ini_simu/00_impt'

path='/lustre/fs22/group/pitz/biaobin/2025/flattop_bcON_1nC/cores_016_impt_32320064_impz_32320064_Np_6.55e+05_phi2_-32.0deg/00_impt/test'


os.chdir(path)

line = 'line'
fname ="lte2.impt"
lte2 = impt(fname,line)

# lte.control["NEMISSION"] = 999

# lte.back2lte(fname="lte3.impt")


#%%
# path = "./ini_simu/01_impz/"
# os.chdir(path)

# line = 'line'
# fname ="lte.impz"
# lte = impz(fname,line)

# # lte.control["NEMISSION"] = 999

# lte.back2lte(fname="lte3.impz")