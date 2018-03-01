#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 11:36:00 2018

@author: FJ
"""

# ----------------------------------------------------------------------------
# Import
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from sympy import plot , symbols , diff, integrate, pi, Rational, Symbol, sin, cos, lambdify, exp, nsimplify, Sum
from sympy.abc import x


# ----------------------------------------------------------------------------
# Parameters
#
N_part=500
dia=1.0

# ----------------------------------------------------------------------------
# Data

# File location
dir = ''
filename='Exercise2_gList_1308.dat'

# Get volume
f = open(dir+filename, "r")
block = f.readlines()
volumes = []
for line in block:
    if line[0:2]=='##':
        words = line.split(sep='\t')
        volumes.append(float(words[-1]))
        position = len(volumes)-1
        print(filename, "heeft op plek ", position, " een V = ", words[-1])
          
data = np.genfromtxt(dir + filename , delimiter ='\t')


# =============================================================================
# rmax=np.sqrt(3)/2*volumes[0]**(1./3.)
# 
# plt.figure()
# 
# Nbins=int(1+3.222*np.log((N_part-1)**2))
# Nbins=45
# 
# # Plot histogram
# plt.hist(data[:,2],Nbins,range=(0,10))
# plt.show()
# 
# # Data histogram
# rho = N_part/volumes[0]
# dr = rmax/Nbins
# 
# counts, bin_edges = np.histogram(data[:,2],bins=np.arange(0,rmax+dr,dr))
# 
# n_id=4*np.pi*N_part/volumes[0]*((bin_edges[:-1]+dr)**3-bin_edges[:-1]**3)/3
# 
# g=(1/N_part)*(counts/n_id)
# 
# gx = np.arange(0,rmax+dr,dr)[:-1]+0.5*dr
# 
# plt.plot(gx,g)
# #plt.xlim(0,5)
# #plt.show()
# =============================================================================

r=data[:,1]
g=data[:,2]

plt.plot(r,g)
plt.show()