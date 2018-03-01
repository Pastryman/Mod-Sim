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
filename='Test.dat'

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

plt.figure()
plt.hist(data[:,2],3)

    
