#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 16:37:20 2019

@author: dori
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
plt.close('all')

rho = 1000.0

ageo = 5.13
bgeo = 0.5
b = 1.0/bgeo#3.0
a = 1.0/ageo**b#rho*np.pi/6.0

def mD(D):
  return a*D**b
def Dm(m):
  return (m/a)**(1/b)

nu=0.0
mu=0.5#1.0/3.0

gam = b*mu
mup = b*nu + b - 1.0

#av = 114.0137
#bv = 0.23437
#
#alpha = 9.292
#beta = 9.623
#ctilde = 622.2
#rho_w = 1000.0
#c = ctilde*(6.0/(np.pi*rho_w))**mu

def lam(N, q):
  return (N*gamma((nu+2)/mu)/(q*gamma((nu+1)/mu)))**mu

def A(N, lam):
  return N*mu*lam**((nu+1)/mu)/gamma((nu+1)/mu)

def LAM(lam):
  return lam*a**mu

def N0(A):
  return A*b*a**(nu+1)

def MkN(N0, LAM, k):
  return N0*gamma((mup+k+1)/gam)/(gam*LAM**((mup+k+1)/gam))

def Mkf(A, lam, k):
  return A*gamma((nu+k+1)/mu)/(mu*lam**((nu+k+1)/mu))

def qN2Dm(q, N, b=3):
  ll = lam(N, q)
  AA = A(N, ll)
  LL = LAM(ll)
  NN = N0(AA)
  return MkN(NN, LL, b+1)/MkN(NN, LL, b)

def LAM2Dm(L):
  return (mup+1.0+b)/(gam*L**(1.0/gam))

def lam2Mm(l):
  #return l**(-1/mu)*(nu+1)/mu
  arg = (nu+1)/mu
  return l**(-1/mu)*gamma(arg+1./mu)/gamma(arg)

def lam2Dm(l):
  return Dm(lam2Mm(l))


q = 1.0e-6
N = np.logspace(0, 7, 1000)

mean_masses = q/N

lams = lam(N, q)
LAMs = LAM(lams)
As = A(N, lams)
M1f_M0f = Mkf(As, lams, 1)/Mkf(As, lams, 0)
Mm = lam2Mm(lams)

plt.figure()
plt.scatter(mean_masses, M1f_M0f)
plt.scatter(mean_masses, Mm)
plt.xlim([mean_masses.min(), mean_masses.max()])
plt.ylim([mean_masses.min(), mean_masses.max()])

D_m = qN2Dm(q, N, b=b)
D_m2 = LAM2Dm(LAMs)
DMm = Dm(mean_masses)
DMm2 = lam2Dm(lams)

plt.figure()
plt.scatter(D_m, DMm)
plt.scatter(D_m, D_m2)
plt.scatter(D_m, DMm2)
plt.xlim([D_m.min(), D_m.max()])
plt.ylim([D_m.min(), D_m.max()])

ratios = DMm/D_m
g1=gamma((nu+2)/mu)
g2=gamma((nu+1)/mu)
g3=gamma((b*nu+2*b+1)/(mu*b))
ratio = (g1/g2)**bgeo*(g1/g3)

