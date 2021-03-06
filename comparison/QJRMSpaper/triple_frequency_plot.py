#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 14:27:31 2019

@author: dori
"""

import matplotlib.pyplot as plt
import netCDF4
import sys
sys.path.append('..')
from READ import slice_data
from READ import read_variables
from statistic import hist_and_plot
import pandas as pd
import numpy as np

rad = pd.read_hdf('../data/radStatCTTmaxDWRregrid.h5', key='stat')
rad.index = pd.to_datetime(rad.index, unit='s')
rpu = pd.read_hdf('../data/radStatCTTmaxDWR.h5', key='stat')
rpu.index = pd.to_datetime(rpu.index, unit='s')
pam = pd.read_hdf('../data/pamStatCTTmaxDWR.h5', key='stat')
pam.index = pd.to_datetime(pam.index, unit='s')

pamtra = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='all_hydro', suffix='pamtra_icon.h5', pamtra=True,
                        varlist=['T', 'Z10', 'Z35', 'Z94', 'unixtime'], minhour=6.0)
radar = read_variables(path='/work/develop/pamtraICON/comparison/data/radar/',
                       hydroset='', suffix='radar_regrid.h5',
                       varlist=['T', 'Z10', 'Z35', 'Z94',
                                'quality_x', 'quality_w',
                                'unixtime'], minhour=6.0)
radar.unixtime = pd.to_datetime(radar.unixtime.astype(np.int64), unit='s')
pamtra.unixtime = pd.to_datetime(pamtra.unixtime.astype(np.int64), unit='s')
rad = rad.loc[radar.unixtime]
rpu = rpu.loc[radar.unixtime]
pam = pam.loc[pamtra.unixtime]

#for col in rpu.columns:
#  radar[col] = rpu[col].values
for col in rad.columns:
  radar[col] = rad[col].values
for col in pam.columns:
  pamtra[col] = pam[col].values
  
data = netCDF4.Dataset('../data/idealized_hydro.nc')
data_simple = netCDF4.Dataset('../data/idealized_hydro_simple.nc')
datavar = data.variables
Zvar = data_simple.variables['Ze']
Ze = Zvar[:]
Ze[:,0,:,2,0,0] = Ze[:,0,:,2,0,0] #+ 10.*np.log10(0.93/0.72)
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

Zc[:,0] -= 0.5
Zr[:,0] -= 0.5

#pamtra = read_prepare(filename='../data/tripex_all_hydro_data_pamtra_icon.h5')
#pamtra[pamtra==9999.0] = np.nan
lognormrule=True
#radar = read_radar(filename='..data/'+campaign+'_data_radar'++'.h5', minhour=minhour)
radarw = slice_data(radar, 'quality_w', maxvalue=8192)
radarx = slice_data(radar, 'quality_x', maxvalue=8192)
radarxw = slice_data(radarw, 'quality_x', maxvalue=8192)

#%% Triple frequency
xlim=[-10, 20]
ylim=[-10, 20]
lw=3
r = hist_and_plot(slice_data(pamtra, 'Z10', minvalue=-2), '3f plot all', yvar='DWRxk', xvar='DWRkw',
              xlabel='DWR Ka W   [dB]', ylabel='DWR X Ka   [dB]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='pamtra3f_all.png', inverty=False, figax=None,
              bins=100,density=True, CFAD=False)
plt.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets', lw=lw)
plt.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals', lw=lw)
plt.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops', lw=lw)
plt.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes', lw=lw)
plt.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel', lw=lw)
plt.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail', lw=lw)
plt.legend(loc=2)
plt.savefig('pamtra3f_all_withcurves.png', dpi=300)

hist_and_plot(slice_data(radarxw, 'Z10', minvalue=5).dropna(subset=['DWRxk', 'DWRkw']),
              '3f plot all', yvar='DWRxk', xvar='DWRkw',
              xlabel='DWR Ka W   [dB]', ylabel='DWR X Ka   [dB]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='radar3f_all.png', inverty=False, figax=None,
              bins=100,density=True, CFAD=False)
plt.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets', lw=lw)
plt.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals', lw=lw)
plt.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops', lw=lw)
plt.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes', lw=lw)
plt.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel', lw=lw)
plt.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail', lw=lw)
plt.savefig('radar3f_all_withcurves.png', dpi=300)

#%% Combined plot for paper
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.5))
ax1.set_aspect('equal')
ax2.set_aspect('equal')
r = hist_and_plot(slice_data(pamtra, 'Z10', minvalue=0), 'Simulated',
                  yvar='DWRxk', xvar='DWRkw', vminmax=(.001, 1),
                  xlabel='DWR$_{K_a W}$   [dB]', ylabel='DWR$_{X K_a}$   [dB]',
                  xlim=xlim, ylim=ylim, lognorm=lognormrule, figax=(f, ax1),
                  savename='pamRad3f_all.png', inverty=False,
                  bins=100, density=True, CFAD=False)

r = hist_and_plot(slice_data(radarxw, 'Z10',
                             minvalue=0).dropna(subset=['DWRxk', 'DWRkw']),
                  'Measured', yvar='DWRxk', xvar='DWRkw', vminmax=(.001, 1),
                  xlabel='DWR$_{K_a W}$   [dB]', ylabel='DWR$_{X K_a}$   [dB]',
                  xlim=xlim, ylim=ylim, lognorm=lognormrule, figax=(f, ax2),
                  savename='pamRad3f_all_no_curves.png', inverty=False,
                  bins=100,#(r[4], r[5]),
                  density=True, CFAD=False)
ax1.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets', lw=lw)
ax1.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals', lw=lw)
ax1.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops', lw=lw)
ax1.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes', lw=lw)
ax1.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel', lw=lw)
ax1.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail', lw=lw)
ax1.legend(loc=2)
ax2.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets', lw=lw)
ax2.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals', lw=lw)
ax2.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops', lw=lw)
ax2.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes', lw=lw)
ax2.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel', lw=lw)
ax2.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail', lw=lw)
f.suptitle('Triple Frequency plots', fontsize=12, fontweight='heavy', y=0.99)
f.tight_layout()
f.savefig('pamRad3f_all.pdf',dpi=600)
f.savefig('pamRad3f_all.png', dpi=300)

#%% Combined plot for paper low maxDWRkw
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.5))
ax1.set_aspect('equal')
ax2.set_aspect('equal')
pam = slice_data(pamtra, 'Z10', minvalue=0.0)
pam = slice_data(pam, 'maxDWRkw', maxvalue=6.0)
r = hist_and_plot(pam, 'Simulated',
                  yvar='DWRxk', xvar='DWRkw', vminmax=(.001, 1),
                  xlabel='DWR$_{K_a W}$   [dB]', ylabel='DWR$_{X K_a}$   [dB]',
                  xlim=xlim, ylim=ylim, lognorm=lognormrule, figax=(f, ax1),
                  savename='dummy.png', inverty=False,
                  bins=100, density=True, CFAD=False)
ax1.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets', lw=lw)
ax1.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals', lw=lw)
ax1.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops', lw=lw)
ax1.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes', lw=lw)
ax1.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel', lw=lw)
ax1.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail', lw=lw)
ax1.legend(loc=2)

rad = slice_data(radarxw, 'Z10', minvalue=0.0).dropna(subset=['DWRxk', 'DWRkw'])
rad = slice_data(rad, 'maxDWRkw', maxvalue=6.0)
r = hist_and_plot(rad, 'Measured',
                  yvar='DWRxk', xvar='DWRkw', vminmax=(.001, 1),
                  xlabel='DWR$_{K_a W}$   [dB]', ylabel='DWR$_{X K_a}$   [dB]',
                  xlim=xlim, ylim=ylim, lognorm=lognormrule, figax=(f, ax2),
                  savename='dummy.png', inverty=False,
                  bins=100,#(r[4], r[5]),
                  density=True, CFAD=False)
ax2.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets', lw=lw)
ax2.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals', lw=lw)
ax2.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops', lw=lw)
ax2.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes', lw=lw)
ax2.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel', lw=lw)
ax2.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail', lw=lw)
f.suptitle('Triple Frequency plots maxDWRkw<6dB', fontsize=12, fontweight='heavy', y=0.99)
f.tight_layout()
f.savefig('pamRad3f_all.png', dpi=300)

#%% Combined plot for paper low maxDWRkw
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.5))
ax1.set_aspect('equal')
ax2.set_aspect('equal')
pam = slice_data(pamtra, 'Z10', minvalue=0.0)
pam = slice_data(pam, 'maxDWRkw', minvalue=6.0)
r = hist_and_plot(pam, 'Simulated',
                  yvar='DWRxk', xvar='DWRkw', vminmax=(.001, 1),
                  xlabel='DWR$_{K_a W}$   [dB]', ylabel='DWR$_{X K_a}$   [dB]',
                  xlim=xlim, ylim=ylim, lognorm=lognormrule, figax=(f, ax1),
                  savename='dummy.png', inverty=False,
                  bins=100, density=True, CFAD=False)
ax1.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets', lw=lw)
ax1.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals', lw=lw)
ax1.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops', lw=lw)
ax1.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes', lw=lw)
ax1.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel', lw=lw)
ax1.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail', lw=lw)
ax1.legend(loc=2)

rad = slice_data(radarxw, 'Z10', minvalue=0.0).dropna(subset=['DWRxk', 'DWRkw'])
rad = slice_data(rad, 'maxDWRkw', minvalue=6.0)
r = hist_and_plot(rad, 'Measured',
                  yvar='DWRxk', xvar='DWRkw', vminmax=(.001, 1),
                  xlabel='DWR$_{K_a W}$   [dB]', ylabel='DWR$_{X K_a}$   [dB]',
                  xlim=xlim, ylim=ylim, lognorm=lognormrule, figax=(f, ax2),
                  savename='dummy.png', inverty=False,
                  bins=100,#(r[4], r[5]),
                  density=True, CFAD=False)
ax2.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets', lw=lw)
ax2.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals', lw=lw)
ax2.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops', lw=lw)
ax2.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes', lw=lw)
ax2.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel', lw=lw)
ax2.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail', lw=lw)
f.suptitle('Triple Frequency plots maxDWRkw>6dB', fontsize=12, fontweight='heavy', y=0.99)
f.tight_layout()
f.savefig('HIGdwrKWpamRad3f_all.png', dpi=300)

#%%

r = hist_and_plot(slice_data(slice_data(pamtra, 'Z10', minvalue=-2), 'T', maxvalue=0),
                      '3f plot T<0', yvar='DWRxk', xvar='DWRkw',
              xlabel='DWR Ka W   [dB]', ylabel='DWR X Ka   [dB]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='pamtra3f_all.png', inverty=False, figax=None,
              bins=100,density=True, CFAD=False)
plt.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets', lw=lw)
plt.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals', lw=lw)
plt.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops', lw=lw)
plt.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes', lw=lw)
plt.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel', lw=lw)
plt.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail', lw=lw)
plt.legend(loc=2)
plt.savefig('pamtra3f_T<0_withcurves.png', dpi=300)

hist_and_plot(slice_data(slice_data(radarxw, 'Z10', minvalue=5), 'T', maxvalue=0).dropna(subset=['DWRxk', 'DWRkw']),
              '3f plot T<0', yvar='DWRxk', xvar='DWRkw',
              xlabel='DWR Ka W   [dB]', ylabel='DWR X Ka   [dB]',
              xlim=xlim, ylim=xlim, lognorm=lognormrule,
              savename='radar3f_all.png', inverty=False, figax=None,
              bins=100,density=True, CFAD=False)
plt.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets', lw=lw)
plt.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals', lw=lw)
plt.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops', lw=lw)
plt.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes', lw=lw)
plt.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel', lw=lw)
plt.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail', lw=lw)
plt.savefig('radar3f_T<0_withcurves.png', dpi=300)

#%%

r = hist_and_plot(slice_data(slice_data(pamtra, 'Z10', minvalue=-2), 'T', minvalue=0),
                      '3f plot T>0', yvar='DWRxk', xvar='DWRkw',
              xlabel='DWR Ka W   [dB]', ylabel='DWR X Ka   [dB]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='pamtra3f_all.png', inverty=False, figax=None,
              bins=100,density=True, CFAD=False)
plt.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets', lw=lw)
plt.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals', lw=lw)
plt.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops', lw=lw)
plt.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes', lw=lw)
plt.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel', lw=lw)
plt.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail', lw=lw)
plt.legend(loc=2)
plt.savefig('pamtra3f_T>0_withcurves.png', dpi=300)

hist_and_plot(slice_data(slice_data(radarxw, 'Z10', minvalue=5), 'T', minvalue=0).dropna(subset=['DWRxk', 'DWRkw']),
              '3f plot T>0', yvar='DWRxk', xvar='DWRkw',
              xlabel='DWR Ka W   [dB]', ylabel='DWR X Ka   [dB]',
              xlim=xlim, ylim=xlim, lognorm=lognormrule,
              savename='radar3f_all.png', inverty=False, figax=None,
              bins=100,density=True, CFAD=False)
plt.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets', lw=lw)
plt.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals', lw=lw)
plt.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops', lw=lw)
plt.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes', lw=lw)
plt.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel', lw=lw)
plt.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail', lw=lw)
plt.savefig('radar3f_T>0_withcurves.png', dpi=300)

#%%

r = hist_and_plot(slice_data(slice_data(pamtra, 'Z10', minvalue=-2), 'T', minvalue=5),
                      '3f plot T>5', yvar='DWRxk', xvar='DWRkw',
              xlabel='DWR Ka W   [dB]', ylabel='DWR X Ka   [dB]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='pamtra3f_all.png', inverty=False, figax=None,
              bins=100,density=True, CFAD=False)
plt.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets', lw=lw)
plt.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals', lw=lw)
plt.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops', lw=lw)
plt.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes', lw=lw)
plt.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel', lw=lw)
plt.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail', lw=lw)
plt.legend(loc=2)
plt.savefig('pamtra3f_T>0_withcurves.png', dpi=300)

hist_and_plot(slice_data(slice_data(radarxw, 'Z10', minvalue=5), 'T', minvalue=5).dropna(subset=['DWRxk', 'DWRkw']),
              '3f plot T>5', yvar='DWRxk', xvar='DWRkw',
              xlabel='DWR Ka W   [dB]', ylabel='DWR X Ka   [dB]',
              xlim=xlim, ylim=xlim, lognorm=lognormrule,
              savename='radar3f_all.png', inverty=False, figax=None,
              bins=100,density=True, CFAD=False)
plt.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets', lw=lw)
plt.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals', lw=lw)
plt.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops', lw=lw)
plt.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes', lw=lw)
plt.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel', lw=lw)
plt.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail', lw=lw)
plt.savefig('radar3f_T>0_withcurves.png', dpi=300)