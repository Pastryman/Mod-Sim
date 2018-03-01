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
# Data (to test)

dir = ''
filename='Exercise2_P50liquid12000.dat'
data = np.genfromtxt(dir + filename , delimiter ='\t')
mean_V=np.mean(data[:,1])
sigma_V=np.std(data[:,1])
print(filename, ": V=",mean_V," ± ",sigma_V)

# ----------------------------------------------------------------------------
# Data (CRYSTAL!)

# Pressure values for CRYSTAL
P=np.array([5,10,15,20,25,30,35,40,45,50])
values=np.zeros((len(P),8))
dir = ''

# Obtaining data from files
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
    
    values[i,0]=P[i]
    values[i,1]=N
    values[i,2]=mean_V
    values[i,3]=sigma_V
    values[i,4]=mean_move_acc
    values[i,5]=sigma_move_acc
    values[i,6]=mean_vol_acc
    values[i,7]=sigma_vol_acc
    
# ----------------------------------------------------------------------------
# Data (LIQUID!)

# Pressure values for LIQUID
P_liq=np.array([5,10,15,20,25,30,35,40,45,50])
values_liq=np.zeros((len(P_liq),8))
dir = ''

# Obtaining data from files
for i in range(0,len(P_liq)):
    filename='Exercise2_P{}_liquid.dat'.format(int(P_liq[i]))
    data = np.genfromtxt(dir + filename , delimiter ='\t')
    
    N=len(data[:,0])
    mean_V=np.mean(data[:,1])
    sigma_V=np.std(data[:,1])
    mean_move_acc=np.mean(data[:,2])
    sigma_move_acc=np.std(data[:,2])
    mean_vol_acc=np.mean(data[:,3])
    sigma_vol_acc=np.std(data[:,3])
    
    values_liq[i,0]=P_liq[i]
    values_liq[i,1]=N
    values_liq[i,2]=mean_V
    values_liq[i,3]=sigma_V
    values_liq[i,4]=mean_move_acc
    values_liq[i,5]=sigma_move_acc
    values_liq[i,6]=mean_vol_acc
    values_liq[i,7]=sigma_vol_acc

# ----------------------------------------------------------------------------
# Plot 1 - Pressure vs Volume

# Plot Data
plt.errorbar(values[:,0], values[:,2],yerr=values[:,4],fmt='.b',label = 'Crystal')
plt.errorbar(values_liq[:,0], values_liq[:,2],yerr=values_liq[:,4],fmt='.r',label = 'Liquid')

# Layout
plt.grid()
plt.xlabel("Dimensionless pressure: $ \\beta P \sigma^{3} $")
plt.ylabel("<Volume> ($\sigma^{3}$)")
plt.ylim([370,700])
plt.legend()

# Visualize
plt.savefig('Fig_Exercise_2.pdf')
plt.show()

# ----------------------------------------------------------------------------
# Plot 2 - Packingfaction vs Volume

# Volume of a particle
vol_part=(4.0/3.0)*np.pi*(dia/2)**3

# Packing fraction - Crystal
pf=N_part*vol_part/values[:,2]
pf_sigma=0.5*abs(N_part*vol_part*(1/(values[:,2]-values[:,3])-1/(values[:,2]+values[:,3])))

# Packing fraction - Liquid
pf_liq=N_part*vol_part/values_liq[:,2]
pf_sigma_liq=0.5*abs(N_part*vol_part*(1/(values_liq[:,2]-values_liq[:,3])-1/(values_liq[:,2]+values_liq[:,3])))

# Errorbar plot for Crystal and Liquid
plt.errorbar(P,pf,yerr=pf_sigma,fmt='.b',label = 'Crystal')
plt.errorbar(P_liq,pf_liq,yerr=pf_sigma_liq,fmt='.r',label = 'Liquid')
plt.grid()

# Including Carnahan and Starling equation of state
# We have to solve the C&S - Equation of State
from sympy.solvers import solve
from sympy import Symbol
eta = Symbol('eta')
bP = Symbol('bP')

# Solve function for eta[P]=...
#f=solve((1+eta+eta**2-eta**3)/(1-eta)**3-bP, eta)              # Remove this '#' to get solutions

# Only use first solution as the others are complex
def etaF(bP):
    f=-(-3*bP + 1)/(3*(bP - 1)) - ((-3*bP + 1)**2/(bP - 1)**2 - 3*(3*bP + 1)/(bP - 1))/(3*((-3*bP + 1)**3/(bP - 1)**3 - 9*(-3*bP + 1)*(3*bP + 1)/(2*(bP - 1)**2) + np.sqrt(-4*((-3*bP + 1)**2/(bP - 1)**2 - 3*(3*bP + 1)/(bP - 1))**3 + (2*(-3*bP + 1)**3/(bP - 1)**3 - 9*(-3*bP + 1)*(3*bP + 1)/(bP - 1)**2 - 27)**2)/2 - 27/2)**(1/3)) - ((-3*bP + 1)**3/(bP - 1)**3 - 9*(-3*bP + 1)*(3*bP + 1)/(2*(bP - 1)**2) + np.sqrt(-4*((-3*bP + 1)**2/(bP - 1)**2 - 3*(3*bP + 1)/(bP - 1))**3 + (2*(-3*bP + 1)**3/(bP - 1)**3 - 9*(-3*bP + 1)*(3*bP + 1)/(bP - 1)**2 - 27)**2)/2 - 27/2)**(1/3)/3
    return f

# Plot Carnahan & Starling EOS
P_ex3=np.linspace(1.8,50,num=1000) # Smaller than packing fraction of 1.8 it gives an error
eta_ex3=etaF(P_ex3)
plt.plot(P_ex3,eta_ex3,"g--",label = 'Carnahan & Starling EOS')

# Layout
plt.xlabel("Dimensionless pressure: $ \\beta P \sigma^{3} $")
plt.ylabel("Packing fraction: $\eta$")
plt.legend()

# Visualize
plt.savefig('Fig_Exercise_3.pdf')
plt.show()

# ----------------------------------------------------------------------------
# Plot 3 - Volume vs Packingfraction (not used in report)

# Formula
def bP(eta):
    f=(1+eta+eta**2-eta**3)/((1-eta)**3)
    return f

# Plot data
plt.errorbar(pf,P,xerr=pf_sigma,fmt='.b',label = 'Crystal')
plt.errorbar(pf_liq,P_liq,xerr=pf_sigma_liq,fmt='.r',label = 'Liquid')

# Plot Formula
eta_ex3=P_ex3=np.linspace(0.35,0.68,num=1000)
plt.plot(eta_ex3,bP(eta_ex3),"g--",label = 'Carnahan & Starling EOS')

# Layout
plt.ylabel("Dimensionless pressure: $ \\beta P \sigma^{3} $")
plt.xlabel("Packing fraction: $\eta$")
plt.legend()
plt.grid()

# Visualize
plt.savefig('Fig_Exercise_3_Inverted.pdf')
plt.show()
