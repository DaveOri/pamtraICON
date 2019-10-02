#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 17:12:11 2019

@author: dori
"""

from scipy.special import gamma
import numpy as np

rho = 1000.0
a = rho*np.pi/6.0
b = 3.0

def mD(D):
  return a*D**b
def Dm(m):
  return (m/a)**(1/b)

nu=0.0
mu=1.0/3.0
gam = b*mu
mup = b*nu + b - 1.0

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

q=0.001
N=1000

Dm1 = Dm(q/N)
ll = lam(N, q)
AA = A(N, ll)
LL = LAM(ll)
NN = N0(AA)
Dm2 = MkN(NN, LL, 4)/MkN(NN, LL, 3)

Dm0 = MkN(NN, LL, 1)/MkN(NN, LL, 0)