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

# Volumes
volumes = [2617,5235,1308,1047,872,747,654,581]


for i in range(0,len(volumes)):
    filename='Exercise2_gList_Volume%d_configs500.dat' % (volumes[i])
          
    data = np.genfromtxt(dir + filename , delimiter ='\t')

    r=data[:,1]
    g=data[:,2]
    
    plt.plot(r,g,'C%d' % (i+1),label='$\\rho=%.3f$' % (500./volumes[i]))
    line=np.ones(len(r))
    plt.plot(r,line,'--k',linewidth=0.5)
    
    plt.xlim(0.8,7)
    plt.ylim(0.5,3.4)
    plt.xlabel("r")
    plt.ylabel("g(r)")
    plt.legend()
    plt.savefig('Fig_gr_v%.0f.pdf' % (volumes[i]))
    plt.show()
    
for i in range(0,len(volumes)):
    filename='Exercise2_gList_Volume%d_configs500.dat' % (volumes[i])
          
    data = np.genfromtxt(dir + filename , delimiter ='\t')

    r=data[:,1]
    g=data[:,2]
    
    plt.plot(r,g,'C%d' % (i+1),label='$\\rho=%.3f$' % (500./volumes[i]))
    
    
line=np.ones(len(r))
plt.plot(r,line,'--k',linewidth=0.5) 
plt.xlim(0.8,7)
plt.ylim(0.5,3.4)
plt.xlabel("r")
plt.ylabel("g(r)")
plt.legend()
plt.savefig('Fig_gr_ALL.pdf')
plt.show()     



