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


volumes=[]

# Measured volumes
volumes = [1308]
filename='Exercise2_gList_Volume%d_configs500.dat' % (volumes[0])
          
data = np.genfromtxt(dir + filename , delimiter ='\t')

r=data[:,1]
g=data[:,2]

plt.plot(r,g,label='V=%.6s' % (volumes[0]))

line=np.ones(len(r))
plt.plot(r,line,'--')
plt.xlim(0,9)
plt.ylim(0.5,1.8)
plt.legend()
plt.show()