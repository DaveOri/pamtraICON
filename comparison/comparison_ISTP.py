#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 16:09:20 2019

@author: dori
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4

from READ import read_prepare
from READ import read_radar
from READ import slice_data
from READ import icon150heights
from statistics import hist_and_plot

data = netCDF4.Dataset('data/idealized_hydro.nc')
data_simple = netCDF4.Dataset('data/idealized_hydro_simple.nc')
datavar = data.variables
Zvar = data_simple.variables['Ze']
Ze = Zvar[:]
Ze[:,0,:,2,0,0] = Ze[:,0,:,2,0,0] + 10.*np.log10(0.93/0.72)
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

campaign = 'tripex'
minhour=6.0

pamtra = read_prepare(hydroset=campaign + '_all_hydro_', minhour=minhour, suffix='_bk')
pamtra[pamtra==9999] = np.nan
lognormrule=True

radar = read_radar(campaign=campaign, minhour=minhour, avg='_avg')
radarw = slice_data(radar, 'quality_w', maxvalue=8192)
radarx = slice_data(radar, 'quality_x', maxvalue=8192)
radarxw = slice_data(radarw, 'quality_x', maxvalue=8192)
#%% Triple frequency
xlim=[-5, 20]
ylim=[-5, 20]
lw=3
hist_and_plot(pamtra, '3f plot all', yvar='DWRxk', xvar='DWRkw',
              xlabel='DWR Ka W   [dB]', ylabel='DWR X Ka   [dB]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='ISTP/pamtra3f_all.png', inverty=False, figax=None,
              bins=100,density=False, CFAD=False)
plt.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets', lw=lw)
plt.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals', lw=lw)
plt.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops', lw=lw)
plt.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes', lw=lw)
plt.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel', lw=lw)
plt.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail', lw=lw)
plt.legend(loc=2)
plt.savefig('ISTP/pamtra3f_all_withcurves.png', dpi=300)

hist_and_plot(radarxw.dropna(subset=['DWRxk', 'DWRkw']),
              '3f plot all', yvar='DWRxk', xvar='DWRkw',
              xlabel='DWR Ka W   [dB]', ylabel='DWR X Ka   [dB]',
              xlim=xlim, ylim=xlim, lognorm=lognormrule,
              savename='ISTP/radar3f_all.png', inverty=False, figax=None,
              bins=100,density=False, CFAD=False)
plt.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets')
plt.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail')
plt.savefig('ISTP/radar3f_all_withcurves.png', dpi=300)

#%% DWR CFAD temperature

hist_and_plot(pamtra, 'CFAD DWRxk', yvar='T', xvar='DWRxk',
              xlabel='DWRxk   [dB]', ylabel='T   [degC]',
              xlim=[-10, 20], ylim=[-20, 0], lognorm=lognormrule,
              savename='ISTP/pamtraCFAD_DWRxk_T.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(pamtra, 'CFAD DWRkw', yvar='T', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='T   [degC]',
              xlim=[-8, 13], ylim=[-20, 0], lognorm=lognormrule,
              savename='ISTP/pamtraCFAD_DWRkw_T.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(radarx.dropna(subset=['DWRxk']),
              'CFAD DWRxk', yvar='T', xvar='DWRxk',
              xlabel='DWRxk   [dB]', ylabel='T   [degC]',
              xlim=[-10, 20], ylim=[-20, 0], lognorm=lognormrule,
              savename='ISTP/radarCFAD_DWRxk_T.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(radarw.dropna(subset=['DWRkw']),
              'CFAD DWRkw', yvar='T', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='T   [degC]',
              xlim=[-8, 13], ylim=[-20, 0], lognorm=lognormrule,
              savename='ISTP/radarCFAD_DWRkw_T.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

#%% Velocity - DWR
inverty = True
xlim = [-5, 15]
ylim = [-3, 0]
lw = 4
h,x,y = hist_and_plot(slice_data(pamtra, 'T', -20, -5),
              'MDV35 vs DWRkw  |  -20<T<-5', yvar='V35', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='MDV Ka [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='ISTP/pamtraSCAT_V_DWRkw_-20<T<-5.png',
              inverty=inverty, figax=None, cmap='jet',
              bins=(100,100),density=False, CFAD=False)
plt.plot(Zi[:,1]-Zi[:,2], -Vi[:,0], label='ice crystals', lw=lw)
plt.plot(Zr[:,1]-Zr[:,2], -Vr[:,0], label='raindrops', lw=lw)
plt.plot(Zs[:,1]-Zs[:,2], -Vs[:,0], label='snowflakes', lw=lw)
plt.plot(Zg[:,1]-Zg[:,2], -Vg[:,0], label='graupel', lw=lw)
plt.plot(Zh[:,1]-Zh[:,2], -Vh[:,0], label='hail', lw=lw)
plt.savefig('ISTP/pamtraSCAT_V_DWRkw_-20<T<-5.png', dpi=300)

Hrad = np.array(sorted(radar['Hgt'].drop_duplicates().values))
Hr_b = np.ndarray((len(Hrad)+1,))
Hr_b[:-1] = Hrad - 15.
Hr_b[-1] = Hrad[-1] + 15. 

h,x,y = hist_and_plot(slice_data(radarw, 'T', -20, -5),
              'MDV35 vs DWRkw  |  -20<T<-5', yvar='V35avg', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='MDV Ka [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='ISTP/radarSCAT_V_DWRkw_-20<T<-5.png',
              inverty=inverty, figax=None, cmap='jet',
              bins=(200,200),density=False, CFAD=False)
plt.plot(Zi[:,1]-Zi[:,2], -Vi[:,0], label='ice crystals', lw=lw)
plt.plot(Zr[:,1]-Zr[:,2], -Vr[:,0], label='raindrops', lw=lw)
plt.plot(Zs[:,1]-Zs[:,2], -Vs[:,0], label='snowflakes', lw=lw)
plt.plot(Zg[:,1]-Zg[:,2], -Vg[:,0], label='graupel', lw=lw)
plt.plot(Zh[:,1]-Zh[:,2], -Vh[:,0], label='hail', lw=lw)
plt.legend(loc=1)
plt.savefig('ISTP/radarSCAT_V_DWRkw_-20<T<-5.png', dpi=300)

#%% Velocity CFAD temperature

h, x, y = hist_and_plot(radar.dropna(subset=['V35avg']), 'MDV Ka',
              yvar='T', xvar='V35avg',
              xlabel='MDV Ka   [m/s]', ylabel='T   [K]',
              xlim=[-5, 2], ylim=[-30, 10], lognorm=False,#vminmax=[0, 6],
              savename='ISTP/radarCFAD_V35_T.png',
              inverty=True, figax=None, CFAD=False,
              bins=(300,300))
plt.savefig('ISTP/radarCFAD_V35_T.png',dpi=300)

hist_and_plot(pamtra, 'MDV Ka',
              yvar='T', xvar='V35',
              xlabel='MDV Ka   [m/s]', ylabel='T   [K]',
              xlim=[-5, 2], ylim=[-30, 10], lognorm=False,# vminmax=[0, 500],
              savename='ISTP/pamtraCFAD_V35_T.png',
              inverty=True, figax=None, CFAD=False,
              bins=(300,300))
plt.savefig('ISTP/pamtraCFAD_V35_T.png',dpi=300)

h, x, y = hist_and_plot(radar.dropna(subset=['V35avg']), 'MDV Ka',
              yvar='T', xvar='V35avg',
              xlabel='MDV Ka   [m/s]', ylabel='T   [K]',
              xlim=[-5, 2], ylim=[-30, 10], lognorm=True,
              savename='ISTP/radarCFAD_V35_T_log.png',
              inverty=True, figax=None, CFAD=False,
              bins=(300,300))
plt.savefig('ISTP/radarCFAD_V35_T_log.png',dpi=300)

hist_and_plot(pamtra, 'MDV Ka',
              yvar='T', xvar='V35',
              xlabel='MDV Ka   [m/s]', ylabel='T   [K]',
              xlim=[-5, 2], ylim=[-30, 10], lognorm=True,
              savename='ISTP/pamtraCFAD_V35_T_log.png',
              inverty=True, figax=None, CFAD=False,
              bins=(600,600))
plt.savefig('ISTP/pamtraCFAD_V35_T_log.png',dpi=300)