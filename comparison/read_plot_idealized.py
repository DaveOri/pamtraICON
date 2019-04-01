#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 10:29:05 2019

@author: dori
"""

import netCDF4
import matplotlib.pyplot as plt

plt.close('all')

data = netCDF4.Dataset('data/idealized_hydro.nc')
data_simple = netCDF4.Dataset('data/idealized_hydro_simple.nc')
datavar = data.variables
Zvar = data_simple.variables['Ze']
Ze = Zvar[:]
V = datavar['Radar_MeanDopplerVel'][:]
S = datavar['Radar_SpectrumWidth'][:]


Zc = Ze[0,0,:,:,0,0]
Zi = Ze[1,0,:,:,0,0]
Zr = Ze[2,0,:,:,0,0]
Zs = Ze[3,0,:,:,0,0]
Zg = Ze[4,0,:,:,0,0]
Zh = Ze[5,0,:,:,0,0]

Vc = V[0,0,:,:,0,0]
Vi = V[1,0,:,:,0,0]
Vr = V[2,0,:,:,0,0]
Vs = V[3,0,:,:,0,0]
Vg = V[4,0,:,:,0,0]
Vh = V[5,0,:,:,0,0]

Sc = S[0,0,:,:,0,0]
Si = S[1,0,:,:,0,0]
Sr = S[2,0,:,:,0,0]
Ss = S[3,0,:,:,0,0]
Sg = S[4,0,:,:,0,0]
Sh = S[5,0,:,:,0,0]


plt.figure()
plt.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets')
plt.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail')
plt.grid()
plt.legend()

plt.figure()
plt.plot(Zc[:,1], Vc[:,1], label='cloud droplets')
plt.plot(Zi[:,1], Vi[:,1], label='ice crystals')
plt.plot(Zr[:,1], Vr[:,1], label='raindrops')
plt.plot(Zs[:,1], Vs[:,1], label='snowflakes')
plt.plot(Zg[:,1], Vg[:,1], label='graupel')
plt.plot(Zh[:,1], Vh[:,1], label='hail')
plt.grid()
plt.legend()

plt.figure()
plt.plot(Zc[:,1], Sc[:,1], label='cloud droplets')
plt.plot(Zi[:,1], Si[:,1], label='ice crystals')
plt.plot(Zr[:,1], Sr[:,1], label='raindrops')
plt.plot(Zs[:,1], Ss[:,1], label='snowflakes')
plt.plot(Zg[:,1], Sg[:,1], label='graupel')
plt.plot(Zh[:,1], Sh[:,1], label='hail')
plt.grid()
plt.legend()

plt.figure()
plt.plot(Zc[:,0]-Zc[:,1], Vc[:,0], label='cloud droplets')
plt.plot(Zi[:,0]-Zi[:,1], Vi[:,0], label='ice crystals')
plt.plot(Zr[:,0]-Zr[:,1], Vr[:,0], label='raindrops')
plt.plot(Zs[:,0]-Zs[:,1], Vs[:,0], label='snowflakes')
plt.plot(Zg[:,0]-Zg[:,1], Vg[:,0], label='graupel')
plt.plot(Zh[:,0]-Zh[:,1], Vh[:,0], label='hail')
plt.grid()
plt.legend()

plt.figure()
plt.plot(Zc[:,1]-Zc[:,2], Vc[:,0], label='cloud droplets')
plt.plot(Zi[:,1]-Zi[:,2], Vi[:,0], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2], Vr[:,0], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2], Vs[:,0], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2], Vg[:,0], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2], Vh[:,0], label='hail')
plt.grid()
plt.legend()

plt.figure()
plt.plot(Zc[:,1]-Zc[:,2], Sc[:,0], label='cloud droplets')
plt.plot(Zi[:,1]-Zi[:,2], Si[:,0], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2], Sr[:,0], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2], Ss[:,0], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2], Sg[:,0], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2], Sh[:,0], label='hail')
plt.grid()
plt.legend()
