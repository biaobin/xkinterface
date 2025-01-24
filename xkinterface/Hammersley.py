# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 03:21:32 2023

@author: lixiangk
"""
from interface import *

class Hammersley:
    def __init__(self, base_in):
        bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                 31, 37, 41, 43, 47, 53, 59, 61, 67,
                 71, 73, 79, 83, 89, 97, 101]
        
        self.base = bases[base_in]
        self.idx = 0

    def set(self, i):
        self.idx = i
        if self.idx < 0:
            self.idx = 0

    def get(self):
        xs = 0
        xsi = 1.0
        i2 = self.idx + 1
        self.idx = i2

        while i2 > 0:
            xsi /= self.base
            i1 = i2 // self.base
            xs += (i2 - self.base * i1) * xsi
            i2 = i1

        return 1-xs

    def __del__(self):
        pass

iprime = 100
ham = Hammersley(25)

fig, ax = plt.subplots(figsize = (5, 4))

#rr = np.array([ham.get() for j in np.arange(1000)])
#ax.hist(rr, bins = 20, histtype = r'step')

#rr = np.array([HaltonNorm2(j, a = -3, b = 0) for j in np.arange(10000)])
rr = np.array([Halton(iprime, j) for j in np.arange(10000)])
ax.hist(rr, bins = 20, histtype = r'step')

rr = HaltonN(10000, iprime = iprime)
ax.hist(rr, bins = 20, histtype = r'step')

