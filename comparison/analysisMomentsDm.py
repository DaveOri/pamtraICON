#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 16:44:40 2019

@author: dori
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import gamma
plt.close('all')

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

av = 114.0137
bv = 0.23437

alpha = 9.292
beta = 9.623
ctilde = 622.2
rho_w = 1000.0
c = ctilde*(6.0/(np.pi*rho_w))**mu

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

def qN2Dm(q, N):
  ll = lam(N, q)
  AA = A(N, ll)
  LL = LAM(ll)
  NN = N0(AA)
  return MkN(NN, LL, 4)/MkN(NN, LL, 3)

def qN2MDVp(q, N):
  ll = lam(N, q)
  AA = A(N, ll)
  return -av*Mkf(AA, ll, 2+bv)/Mkf(AA, ll, 2)


def qN2SWp(q, N):
  ll = lam(N, q)
  AA = A(N, ll)
  mdv = -qN2MDVp(q, N)
  return np.sqrt(av*av*Mkf(AA, ll, 2+2*bv)/Mkf(AA, ll, 2) - mdv*mdv)

def qN2MDVa(q, N):
  ll = lam(N, q)
  AA = A(N, ll)
  return -(alpha - beta*Mkf(AA, ll+c, 2)/Mkf(AA, ll, 2))

def qN2SWa(q, N):
  ll = lam(N, q)
  AA = A(N, ll)
  mdv = -qN2MDVa(q, N)
  return np.sqrt(alpha*alpha + beta*(-2*alpha*Mkf(AA, ll+c, 2) + beta*Mkf(AA, ll+2*c, 2))/Mkf(AA, ll, 2) - mdv*mdv)

def qN2SKa(q, N):
  ll = lam(N, q)
  AA = A(N, ll)
  mdv = -qN2MDVa(q, N)
  sw = qN2SWa(q, N)
  M3 = alpha**3*Mkf(AA, ll, 2) - 3*alpha**2*beta*Mkf(AA, ll+c, 2) + \
       3*alpha*beta**2*Mkf(AA, ll+2*c, 2) - beta**3*Mkf(AA, ll+3*c, 2)
  return M3/(sw**3*Mkf(AA, ll, 2)) - 3*mdv/sw - (mdv/sw)**3
  
def qN2moments(q, N):
  ll = lam(N, q)
  AA = 1#A(N, ll)
  M2 = Mkf(AA, ll, 2)
  M2_b = Mkf(AA, ll, 2+bv)
  mdv = av*M2_b/M2
  M2b_2 = Mkf(AA, ll, 2*bv+2)
  sw = np.sqrt(av**2*M2b_2/M2 - mdv*mdv)
  M2_3b = Mkf(AA, ll, 2+3*bv)
  sk = av**3*M2_3b/(sw**3*Mkf(AA, ll, 2)) - 3*mdv/sw - (mdv/sw)**3
  return -mdv, sw, sk

q = 1.0e-6
N = np.logspace(0, 7, 1000)

mdvp, swp, skp = qN2moments(q, N)
print(skp[0])
mdva = qN2MDVa(q, N)
swa = qN2SWa(q, N)
ska = qN2SKa(q, N)

D = qN2Dm(q, N)*1000

plt.figure()
plt.plot(D, mdvp)
plt.plot(D, mdva)

plt.figure()
plt.plot(D, swp)
plt.plot(D, swa)

plt.figure()
plt.plot(D, skp)
plt.plot(D, ska)