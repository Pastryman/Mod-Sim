#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 11:36:00 2018

@author: FJ
"""

# -----------------
# Import
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from sympy import plot , symbols , diff, integrate, pi, Rational, Symbol, sin, cos, lambdify, exp, nsimplify, Sum
from sympy.abc import x


# -----------------
# Parameters
#
N_part=500
dia=1.0


# -----------------
# Data
#

P=np.array([15,20,25,30,35,40,45,50])
waarden=np.zeros((len(P),8))
dir = '/Users/FJ/C/GitHub/Mod-Sim/Week 3/'

for i in range(0,len(P)):
    filename='Exercise2_P{}.dat'.format(int(P[i]))
    data = np.genfromtxt(dir + filename , delimiter ='\t')
    
    N=len(data[:,0])
    mean_V=np.mean(data[:,1])
    sigma_V=np.std(data[:,1])
    mean_move_acc=np.mean(data[:,2])
    sigma_move_acc=np.std(data[:,2])
    mean_vol_acc=np.mean(data[:,3])
    sigma_vol_acc=np.std(data[:,3])
    
    waarden[i,0]=P[i]
    waarden[i,1]=N
    waarden[i,2]=mean_V
    waarden[i,3]=sigma_V
    waarden[i,4]=mean_move_acc
    waarden[i,5]=sigma_move_acc
    waarden[i,6]=mean_vol_acc
    waarden[i,7]=sigma_vol_acc

# -----------------
# Beeldgeven
#


# Pressure vs Volume
plt.errorbar(waarden[:,0], waarden[:,2],yerr=waarden[:,4],fmt='.')
plt.grid()
plt.xlabel("Pressure ($ \\beta P \sigma^{3} $)")
plt.ylabel("Volume ($\sigma^{3}$)")
plt.ylim([380,500])
#plt.title("Loglog-plot of Mean Square Deviation of $N_{hit}/N$")
plt.show()

# Packingfaction vs Volume
vol_part=(4.0/3.0)*np.pi*(dia/2)**3
pf=N_part*vol_part/waarden[:,2]
pf_sigma=0.5*abs(N_part*vol_part*(1/(waarden[:,2]-waarden[:,3])-1/(waarden[:,2]+waarden[:,3])))

plt.errorbar(P,pf,yerr=pf_sigma,fmt='.r')
plt.grid()
#plt.xlabel("Pressure (beta P sigma^{3})")
#plt.ylabel("Packing fraction (eta)")

plt.show()

#plt.savefig('Fig - Exercise 1.pdf')



# -----------------
# Exercise 3
#

eta=np.linspace(0.5,0.75,num=1000)
def formule(eta):
    f=(1+eta+eta**2-eta**3)/(1-eta)**3
    return f

plt.plot(eta,formule(eta))
plt.plot(pf,formule(pf),".r")
plt.grid()
#plt.yscale('log')
plt.show()
