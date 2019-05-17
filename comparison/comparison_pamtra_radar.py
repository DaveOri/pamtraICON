#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 09:28:57 2019

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

pamtra = read_prepare(hydroset=campaign + '_all_hydro_', minhour=minhour)
pamtraTl0 = slice_data(pamtra, 'T', maxvalue=0)
lognormrule=True

#%% CFAD Sensitivity
xlim = [-60, 50]
ylim = [0, 12000]
h1,x1,y1 = hist_and_plot(pamtra, 'sensitivity cone X', yvar='Hgt', xvar='Z10',
              xlabel='Zx   [dBZ]', ylabel='Hgt   [m]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraCFAD_Zx_Hgt.png',
              inverty=False, figax=None,
              bins=(100,icon150heights[::-1]),
              density=False, CFAD=False)
#p = pamtra[['Z10', 'Hgt']].dropna()
#hst, xedge, yedge = np.histogram2d(p['Z10'], p['Hgt'], bins=(100,icon150heights[::-1]))

hist_and_plot(pamtra, 'sensitivity cone Ka', yvar='Hgt', xvar='Z35',
              xlabel='Zk   [dBZ]', ylabel='Hgt   [m]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraCFAD_Zk_Hgt.png',
              inverty=False, figax=None,
              bins=(100,icon150heights[::-1]))

hist_and_plot(pamtra, 'sensitivity cone W', yvar='Hgt', xvar='Z94',
              xlabel='Zx   [dBZ]', ylabel='Hgt   [m]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraCFAD_Zw_Hgt.png',
              inverty=False, figax=None,
              bins=(100,icon150heights[::-1]))

hist_and_plot(pamtra, 'MDV Ka',
              yvar='T', xvar='V35',
              xlabel='MDV Ka   [m/s]', ylabel='T   [K]',
              xlim=[-5, 2], ylim=[-30, 10], lognorm=lognormrule,
              savename='tripex/3f/pamtraCFAD_V35_T.png',
              inverty=True, figax=None,
              bins=(100,150))

hist_and_plot(pamtra, 'MDV X',
              yvar='T', xvar='V10',
              xlabel='MDV X   [m/s]', ylabel='T   [K]',
              xlim=[-5, 2], ylim=[-30, 10], lognorm=lognormrule,
              savename='tripex/3f/pamtraCFAD_V10_T.png',
              inverty=True, figax=None,
              bins=(100,150))

hist_and_plot(pamtra, 'MDV W',
              yvar='T', xvar='V94',
              xlabel='MDV W   [m/s]', ylabel='T   [K]',
              xlim=[-5, 2], ylim=[-30, 10], lognorm=lognormrule,
              savename='tripex/3f/pamtraCFAD_V94_T.png',
              inverty=True, figax=None,
              bins=(100,150))

radar = read_radar(campaign=campaign, minhour=minhour)
Hrad = np.array(sorted(radar['Hgt'].drop_duplicates().values))
Hr_b = np.ndarray((len(Hrad)+1,))
Hr_b[:-1] = Hrad - 15.
Hr_b[-1] = Hrad[-1] + 15. 

hist_and_plot(radar.dropna(subset=['Z35']), 'sensitivity cone Ka',
              yvar='Hgt', xvar='Z35',
              xlabel='Zk   [dBZ]', ylabel='Hgt   [m]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarCFAD_Zk_Hgt.png',
              inverty=False, figax=None,
              bins=(100,Hr_b))

hist_and_plot(radar.dropna(subset=['V35']), 'MDV Ka',
              yvar='T', xvar='V35',
              xlabel='MDV Ka   [m/s]', ylabel='T   [K]',
              xlim=[-5, 2], ylim=[-30, 10], lognorm=lognormrule,
              savename='tripex/3f/radarCFAD_V35_T.png',
              inverty=True, figax=None,
              bins=(300,300))

hist_and_plot(radar.dropna(subset=['V10']), 'MDV X',
              yvar='T', xvar='V10',
              xlabel='MDV X   [m/s]', ylabel='T   [K]',
              xlim=[-5, 2], ylim=[-30, 10], lognorm=lognormrule,
              savename='tripex/3f/radarCFAD_V10_T.png',
              inverty=True, figax=None,
              bins=(150,200))

hist_and_plot(radar.dropna(subset=['V94']), 'MDV W',
              yvar='T', xvar='V94',
              xlabel='MDV W   [m/s]', ylabel='T   [K]',
              xlim=[-5, 2], ylim=[-30, 10], lognorm=lognormrule,
              savename='tripex/3f/radarCFAD_V94_T.png',
              inverty=True, figax=None,
              bins=(800,200))

radar = slice_data(radar, 'quality_x', maxvalue=8192)
hist_and_plot(radar.dropna(subset=['Z10']), 'sensitivity cone X',
              yvar='Hgt', xvar='Z10',
              xlabel='Zx   [dBZ]', ylabel='Hgt   [m]',
              xlim=[-60, 50], ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarCFAD_Zx_Hgt.png',
              inverty=False, figax=None,
              bins=(100,Hr_b), density=False, CFAD=False)


radar = read_radar(campaign=campaign, minhour=minhour)
radar = slice_data(radar, 'quality_w', maxvalue=8192)
hist_and_plot(radar.dropna(subset=['Z94']), 'sensitivity cone W',
              yvar='Hgt', xvar='Z94',
              xlabel='Zx   [dBZ]', ylabel='Hgt   [m]',
              xlim=[-60, 50], ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarCFAD_Zw_Hgt.png',
              inverty=False, figax=None,
              bins=(100,Hr_b))
radar = slice_data(radar, 'quality_x', maxvalue=8192)
radarTl0 = slice_data(radar, 'T', maxvalue=0)

plt.close('all')

#%% Triple frequency

hist_and_plot(pamtra, '3f plot all', yvar='DWRxk', xvar='DWRkw',
              xlabel='DWR Ka W   [dB]', ylabel='DWR X Ka   [dB]',
              xlim=[-5, 20], ylim=[-5, 20], lognorm=lognormrule,
              savename='tripex/3f/pamtra3f_all.png', inverty=False, figax=None,
              bins=100,density=True, CFAD=True)
plt.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets')
plt.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail')
plt.legend()
plt.savefig('tripex/pamtra3f_all_withcurves.png')

hist_and_plot(slice_data(pamtra, 'T', maxvalue=0),
              '3f plot T<0', yvar='DWRxk', xvar='DWRkw',
              xlabel='DWR Ka W   [dB]', ylabel='DWR X Ka   [dB]',
              xlim=[-5, 20], ylim=[-5, 20], lognorm=lognormrule,
              savename='tripex/3f/pamtra3f_T<0.png', inverty=False, figax=None,
              bins=100,density=True, CFAD=True)

hist_and_plot(radar.dropna(subset=['DWRxk', 'DWRkw']),
              '3f plot all', yvar='DWRxk', xvar='DWRkw',
              xlabel='DWR Ka W   [dB]', ylabel='DWR X Ka   [dB]',
              xlim=[-5, 20], ylim=[-5, 20], lognorm=lognormrule,
              savename='tripex/3f/radar3f_all.png', inverty=False, figax=None,
              bins=100,density=True, CFAD=True)
plt.plot(Zc[:,1]-Zc[:,2],Zc[:,0]-Zc[:,1], label='cloud droplets')
plt.plot(Zi[:,1]-Zi[:,2],Zi[:,0]-Zi[:,1], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2],Zr[:,0]-Zr[:,1], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2],Zs[:,0]-Zs[:,1], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2],Zg[:,0]-Zg[:,1], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2],Zh[:,0]-Zh[:,1], label='hail')
plt.savefig('tripex/radar3f_all_withcurves.png')

hist_and_plot(slice_data(radar.dropna(subset=['DWRxk', 'DWRkw']), 'T', maxvalue=0),
              '3f plot T<0', yvar='DWRxk', xvar='DWRkw',
              xlabel='DWR Ka W   [dB]', ylabel='DWR X Ka   [dB]',
              xlim=[-5, 20], ylim=[-5, 20], lognorm=lognormrule,
              savename='tripex/3f/radar3f_T<0.png', inverty=False, figax=None,
              bins=100,density=True, CFAD=True)

#%% DWR CFAD temperature

hist_and_plot(pamtra, 'CFAD DWRxk', yvar='T', xvar='DWRxk',
              xlabel='DWRxk   [dB]', ylabel='T   [degC]',
              xlim=[-10, 20], ylim=[-20, 0], lognorm=lognormrule,
              savename='tripex/3f/pamtraCFAD_DWRxk_T.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(pamtra, 'CFAD DWRkw', yvar='T', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='T   [degC]',
              xlim=[-8, 13], ylim=[-20, 0], lognorm=lognormrule,
              savename='tripex/3f/pamtraCFAD_DWRkw_T.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

radar = read_radar(campaign=campaign, minhour=minhour)
radar = slice_data(radar, 'quality_x', maxvalue=8192)
hist_and_plot(radar.dropna(subset=['DWRxk']),
              'CFAD DWRxk', yvar='T', xvar='DWRxk',
              xlabel='DWRxk   [dB]', ylabel='T   [degC]',
              xlim=[-10, 20], ylim=[-20, 0], lognorm=lognormrule,
              savename='tripex/3f/radarCFAD_DWRxk_T.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

radar = read_radar(campaign=campaign, minhour=minhour)
radar = slice_data(radar, 'quality_w', maxvalue=8192)
hist_and_plot(radar.dropna(subset=['DWRkw']),
              'CFAD DWRkw', yvar='T', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='T   [degC]',
              xlim=[-8, 13], ylim=[-20, 0], lognorm=lognormrule,
              savename='tripex/3f/radarCFAD_DWRkw_T.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

plt.close('all')

#%% V-Z plots for each DWR range

hist_and_plot(slice_data(pamtra, 'DWRkw', 0, 2),
              'V10 vs Z10  |  0<DWRkw<2', yvar='V10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='MDV X [m/s]',
              xlim=[-60, 50], ylim=[-10, 0], lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_V_Z_DWRkw_0_2.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(slice_data(pamtra, 'DWRkw', 2, 4),
              'V10 vs Z10  |  2<DWRkw<4', yvar='V10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='MDV X [m/s]',
              xlim=[-60, 50], ylim=[-10, 0], lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_V_Z_DWRkw_2_4.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(slice_data(pamtra, 'DWRxk', 4, 6),
              'V10 vs Z10  |  4<DWRkw<6', yvar='V10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='MDV X [m/s]',
              xlim=[-60, 50], ylim=[-10, 0], lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_V_Z_DWRkw_4_6.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(slice_data(pamtra, 'DWRkw', 6, 8),
              'V10 vs Z10  |  6<DWRxk<8', yvar='V10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='MDV X [m/s]',
              xlim=[-60, 50], ylim=[-10, 0], lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_V_Z_DWRkw_6_8.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)


hist_and_plot(slice_data(pamtraTl0, 'DWRkw', 0, 2),
              'V10 vs Z10  |  0<DWRkw<2  |  T<0', yvar='V10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='MDV X [m/s]',
              xlim=[-60, 50], ylim=[-10, 0], lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_V_Z_DWRkw_0_2_Tneg.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(slice_data(pamtraTl0, 'DWRkw', 2, 4),
              'V10 vs Z10  |  2<DWRkw<4  |  T<0', yvar='V10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='MDV X [m/s]',
              xlim=[-60, 50], ylim=[-10, 0], lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_V_Z_DWRkw_2_4_Tneg.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

radar = read_radar(campaign=campaign, minhour=minhour)
radar = slice_data(radar, 'quality_x', maxvalue=8192)
hist_and_plot(slice_data(radar.dropna(subset=['Z10','V10']), 'DWRkw', 0, 2),
              'V10 vs Z10  |  0<DWRkw<2', yvar='V10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='MDV X [m/s]',
              xlim=[-60, 50], ylim=[-10, 0], lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_V_Z_DWRkw_0_2.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(slice_data(radar.dropna(subset=['Z10','V10']), 'DWRkw', 2, 4),
              'V10 vs Z10  |  2<DWRkw<4', yvar='V10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='MDV X [m/s]',
              xlim=[-60, 50], ylim=[-10, 0], lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_V_Z_DWRkw_2_4.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(slice_data(radar.dropna(subset=['Z10','V10']), 'DWRxk', 4, 6),
              'V10 vs Z10  |  4<DWRkw<6', yvar='V10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='MDV X [m/s]',
              xlim=[-60, 50], ylim=[-10, 0], lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_V_Z_DWRkw_4_6.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(slice_data(radar.dropna(subset=['Z10','V10']), 'DWRkw', 6, 8),
              'V10 vs Z10  |  6<DWRxk<8', yvar='V10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='MDV X [m/s]',
              xlim=[-60, 50], ylim=[-10, 0], lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_V_Z_DWRkw_6_8.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)


hist_and_plot(slice_data(radarTl0.dropna(subset=['Z10','V10']), 'DWRkw', 0, 2),
              'V10 vs Z10  |  0<DWRkw<2  |  T<0', yvar='V10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='MDV X [m/s]',
              xlim=[-60, 50], ylim=[-10, 0], lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_V_Z_DWRkw_0_2_Tneg.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(slice_data(radarTl0.dropna(subset=['Z10','V10']), 'DWRkw', 2, 4),
              'V10 vs Z10  |  2<DWRkw<4  |  T<0', yvar='V10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='MDV X [m/s]',
              xlim=[-60, 50], ylim=[-10, 0], lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_V_Z_DWRkw_2_4_Tneg.png',
              inverty=True, figax=None,
              bins=(100,100),density=True, CFAD=True)

plt.close('all')

#%% SW-Z plots for each DWR range

inverty = False
xlim = [-60, 50]
ylim = [0, 2]

hist_and_plot(slice_data(pamtra, 'DWRkw', 0, 2),
              'SW10 vs Z10  |  0<DWRkw<2', yvar='W10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_W_Z_DWRkw_0_2.png',
              inverty=inverty, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(slice_data(pamtra, 'DWRkw', 2, 4),
              'SW10 vs Z10  |  2<DWRkw<4', yvar='W10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_W_Z_DWRkw_2_4.png',
              inverty=inverty, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(slice_data(pamtra, 'DWRxk', 4, 6),
              'SW10 vs Z10  |  4<DWRkw<6', yvar='W10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_W_Z_DWRkw_4_6.png',
              inverty=inverty, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(slice_data(pamtra, 'DWRkw', 6, 8),
              'SW10 vs Z10  |  6<DWRxk<8', yvar='W10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_W_Z_DWRkw_6_8.png',
              inverty=inverty, figax=None,
              bins=(100,100),density=True, CFAD=True)


hist_and_plot(slice_data(pamtraTl0, 'DWRkw', 0, 2),
              'SW10 vs Z10  |  0<DWRkw<2  |  T<0', yvar='W10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_W_Z_DWRkw_0_2_Tneg.png',
              inverty=inverty, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(slice_data(pamtraTl0, 'DWRkw', 2, 4),
              'SW10 vs Z10  |  2<DWRkw<4  |  T<0', yvar='W10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_W_Z_DWRkw_2_4_Tneg.png',
              inverty=inverty, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(slice_data(radar.dropna(subset=['Z10','W10']), 'DWRkw', 0, 2),
              'SW10 vs Z10  |  0<DWRkw<2', yvar='W10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_W_Z_DWRkw_0_2.png',
              inverty=inverty, figax=None,
              bins=(65,65),density=True, CFAD=True)

hist_and_plot(slice_data(radar.dropna(subset=['Z10','W10']), 'DWRkw', 2, 4),
              'SW10 vs Z10  |  2<DWRkw<4', yvar='W10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_W_Z_DWRkw_2_4.png',
              inverty=inverty, figax=None,
              bins=(60,60),density=True, CFAD=True)

hist_and_plot(slice_data(radar.dropna(subset=['Z10','W10']), 'DWRxk', 4, 6),
              'SW10 vs Z10  |  4<DWRkw<6', yvar='W10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_W_Z_DWRkw_4_6.png',
              inverty=inverty, figax=None,
              bins=(60,60),density=True, CFAD=True)

hist_and_plot(slice_data(radar.dropna(subset=['Z10','W10']), 'DWRkw', 6, 8),
              'SW10 vs Z10  |  6<DWRxk<8', yvar='W10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_W_Z_DWRkw_6_8.png',
              inverty=inverty, figax=None,
              bins=(60,60),density=True, CFAD=True)


hist_and_plot(slice_data(radarTl0.dropna(subset=['Z10','W10']), 'DWRkw', 0, 2),
              'SW10 vs Z10  |  0<DWRkw<2  |  T<0', yvar='W10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_W_Z_DWRkw_0_2_Tneg.png',
              inverty=inverty, figax=None,
              bins=(45,45),density=True, CFAD=True)

hist_and_plot(slice_data(radarTl0.dropna(subset=['Z10','W10']), 'DWRkw', 2, 4),
              'SW10 vs Z10  |  2<DWRkw<4  |  T<0', yvar='W10', xvar='Z10',
              xlabel='Z X   [dBZ]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_W_Z_DWRkw_2_4_Tneg.png',
              inverty=inverty, figax=None,
              bins=(50,50),density=True, CFAD=True)

plt.close('all')

#%% V-DWR plots for each T range

inverty = True
xlim = [-5, 15]
ylim = [-5, 0]

radar = read_radar(campaign=campaign, minhour=minhour)
radar = slice_data(radar, 'quality_x', maxvalue=8192)
radar = slice_data(radar, 'quality_w', maxvalue=8192)
h,x,y = hist_and_plot(slice_data(pamtra, 'T', 0),
              'V10 vs DWRkw  |  T>0', yvar='V10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='V X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_V_DWRkw_T>0.png',
              inverty=inverty, figax=None,
              bins=(100,100),density=True, CFAD=True)
plt.plot(Zi[:,1]-Zi[:,2], -Vi[:,0], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2], -Vr[:,0], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2], -Vs[:,0], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2], -Vg[:,0], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2], -Vh[:,0], label='hail')
#plt.plot(Zc[:,1]-Zc[:,2], -Vc[:,0], label='cloud droplets')

plt.legend()
plt.savefig('tripex/pamtraSCAT_V_DWRkw_T>0withcurves.png')

hist_and_plot(slice_data(radar, 'T', 0),
              'V10 vs DWRkw  |  T>0', yvar='V10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='V X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_V_DWRkw_T>0.png',
              inverty=inverty, figax=None,
              bins=(x,y),density=True, CFAD=True)
plt.plot(Zi[:,1]-Zi[:,2], -Vi[:,0], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2], -Vr[:,0], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2], -Vs[:,0], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2], -Vg[:,0], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2], -Vh[:,0], label='hail')
#plt.plot(Zc[:,1]-Zc[:,2], -Vc[:,0], label='cloud droplets')

#plt.legend()
plt.savefig('tripex/radarSCAT_V_DWRkw_T>0withcurves.png')

radar = read_radar(campaign=campaign, minhour=minhour)
radar = slice_data(radar, 'quality_x', maxvalue=8192)
radar = slice_data(radar, 'quality_w', maxvalue=8192)
h,x,y = hist_and_plot(slice_data(pamtra, 'T', -5, 0),
              'V10 vs DWRkw  |  -5<T<0', yvar='V10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='V X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_V_DWRkw_T_-5_0.png',
              inverty=inverty, figax=None,
              bins=(100,100),density=True, CFAD=True)
plt.plot(Zi[:,1]-Zi[:,2], -Vi[:,0], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2], -Vr[:,0], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2], -Vs[:,0], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2], -Vg[:,0], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2], -Vh[:,0], label='hail')
#plt.plot(Zc[:,1]-Zc[:,2], -Vc[:,0], label='cloud droplets')

plt.legend()
plt.savefig('tripex/pamtraSCAT_V_DWRkw_-5_0withcurves.png')

hist_and_plot(slice_data(radar, 'T', -5, 0),
              'V10 vs DWRkw  |  -5<T<0', yvar='V10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='V X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_V_DWRkw_T_-5_0.png',
              inverty=inverty, figax=None,
              bins=(x,y),density=True, CFAD=True)
plt.plot(Zi[:,1]-Zi[:,2], -Vi[:,0], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2], -Vr[:,0], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2], -Vs[:,0], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2], -Vg[:,0], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2], -Vh[:,0], label='hail')
#plt.plot(Zc[:,1]-Zc[:,2], -Vc[:,0], label='cloud droplets')

#plt.legend()
plt.savefig('tripex/radarSCAT_V_DWRkw_-5_0withcurves.png')

h,x,y = hist_and_plot(slice_data(pamtra, 'T', -10, -5),
              'V10 vs DWRkw  |  -10<T<-5', yvar='V10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='V X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_V_DWRkw_T_-10_-5.png',
              inverty=inverty, figax=None,
              bins=(100,100),density=True, CFAD=True)
plt.plot(Zi[:,1]-Zi[:,2], -Vi[:,0], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2], -Vr[:,0], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2], -Vs[:,0], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2], -Vg[:,0], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2], -Vh[:,0], label='hail')
plt.legend()
plt.savefig('tripex/pamtraSCAT_V_DWRkw_-10_-5withcurves.png')

hist_and_plot(slice_data(radar.dropna(), 'T', -10, -5),
              'V10 vs DWRkw  |  -10<T<-5', yvar='V10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='V X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_V_DWRkw_T_-10_-5.png',
              inverty=inverty, figax=None,
              bins=(80,80),density=True, CFAD=True)
plt.plot(Zi[:,1]-Zi[:,2], -Vi[:,0], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2], -Vr[:,0], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2], -Vs[:,0], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2], -Vg[:,0], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2], -Vh[:,0], label='hail')
plt.savefig('tripex/radarSCAT_V_DWRkw_-10_-5withcurves.png')

h,x,y = hist_and_plot(slice_data(pamtra, 'T', -15, -10),
              'V10 vs DWRkw  |  -15<T<-10', yvar='V10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='V X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_V_DWRkw_T_-15_-10.png',
              inverty=inverty, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(slice_data(radar.dropna(), 'T', -15, -10),
              'V10 vs DWRkw  |  -15<T<-10', yvar='V10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='V X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_V_DWRkw_T_-15_-10.png',
              inverty=inverty, figax=None,
              bins=(72,72),density=True, CFAD=True)

h,x,y = hist_and_plot(slice_data(pamtra.dropna(), 'T', -20, -15),
              'V10 vs DWRkw  |  -20<T<-15', yvar='V10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='V X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_V_DWRkw_T_-20_-15.png',
              inverty=inverty, figax=None,
              bins=(100,100),density=True, CFAD=True)

hist_and_plot(slice_data(radar.dropna(), 'T', -20, -15),
              'V10 vs DWRkw  |  -20<T<-15', yvar='V10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='V X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_V_DWRkw_T_-20_-15.png',
              inverty=inverty, figax=None,
              bins=(40,40),density=True, CFAD=True)

plt.close('all')

#%% V-DWR plots for each T range

inverty = False
xlim = [-5, 15]
ylim = [0, 2]

radar = read_radar(campaign=campaign, minhour=minhour)
radar = slice_data(radar, 'quality_x', maxvalue=8192)
radar = slice_data(radar, 'quality_w', maxvalue=8192)

h,x,y = hist_and_plot(slice_data(pamtra, 'T', 0),
              'SW10 vs DWRkw  |  T>0', yvar='W10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_SW_DWRkw_T>0.png',
              inverty=inverty, figax=None,
              bins=(80,80),density=True, CFAD=True)

hist_and_plot(slice_data(radar, 'T', 0),
              'SW10 vs DWRkw  |  T>0', yvar='W10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_SW_DWRkw_T>0.png',
              inverty=inverty, figax=None,
              bins=(x,y),density=True, CFAD=True)

radar = read_radar(campaign=campaign, minhour=minhour)
radar = slice_data(radar, 'quality_x', maxvalue=8192)
radar = slice_data(radar, 'quality_w', maxvalue=8192)
h,x,y = hist_and_plot(slice_data(pamtra, 'T', -5, 0),
              'SW10 vs DWRkw  |  -5<T<0', yvar='W10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_SW_DWRkw_T_-5_0.png',
              inverty=inverty, figax=None,
              bins=(80,80),density=True, CFAD=True)
plt.plot(Zi[:,1]-Zi[:,2], Si[:,0], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2], Sr[:,0], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2], Ss[:,0], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2], Sg[:,0], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2], Sh[:,0], label='hail')
plt.legend()
plt.savefig('tripex/pamtraSCAT_SW_DWRkw_-5_0withcurves.png')

hist_and_plot(slice_data(radar, 'T', -5, 0),
              'SW10 vs DWRkw  |  -5<T<0', yvar='W10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_SW_DWRkw_T_-5_0.png',
              inverty=inverty, figax=None,
              bins=(x,y),density=True, CFAD=True)
plt.plot(Zi[:,1]-Zi[:,2], Si[:,0], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2], Sr[:,0], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2], Ss[:,0], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2], Sg[:,0], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2], Sh[:,0], label='hail')
plt.savefig('tripex/radarSCAT_SW_DWRkw_-5_0withcurves.png')

h,x,y = hist_and_plot(slice_data(pamtra, 'T', -10, -5),
              'SW10 vs DWRkw  |  -10<T<-5', yvar='W10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_SW_DWRkw_T_-10_-5.png',
              inverty=inverty, figax=None,
              bins=(80,80),density=True, CFAD=True)
plt.plot(Zi[:,1]-Zi[:,2], Si[:,0], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2], Sr[:,0], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2], Ss[:,0], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2], Sg[:,0], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2], Sh[:,0], label='hail')
plt.legend()
plt.savefig('tripex/pamtraSCAT_SW_DWRkw_-10_-5withcurves.png')

hist_and_plot(slice_data(radar.dropna(), 'T', -10, -5),
              'SW10 vs DWRkw  |  -10<T<-5', yvar='W10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_SW_DWRkw_T_-10_-5.png',
              inverty=inverty, figax=None,
              bins=(x,y),density=True, CFAD=True)
plt.plot(Zi[:,1]-Zi[:,2], Si[:,0], label='ice crystals')
plt.plot(Zr[:,1]-Zr[:,2], Sr[:,0], label='raindrops')
plt.plot(Zs[:,1]-Zs[:,2], Ss[:,0], label='snowflakes')
plt.plot(Zg[:,1]-Zg[:,2], Sg[:,0], label='graupel')
plt.plot(Zh[:,1]-Zh[:,2], Sh[:,0], label='hail')
plt.savefig('tripex/radarSCAT_SW_DWRkw_-10_-5withcurves.png')

h,x,y = hist_and_plot(slice_data(pamtra, 'T', -15, -10),
              'SW10 vs DWRkw  |  -15<T<-10', yvar='W10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/pamtraSCAT_SW_DWRkw_T_-15_-10.png',
              inverty=inverty, figax=None,
              bins=(80,80),density=True, CFAD=True)

hist_and_plot(slice_data(radar.dropna(), 'T', -15, -10),
              'SW10 vs DWRkw  |  -15<T<-10', yvar='W10', xvar='DWRkw',
              xlabel='DWRkw   [dB]', ylabel='SW X [m/s]',
              xlim=xlim, ylim=ylim, lognorm=lognormrule,
              savename='tripex/3f/radarSCAT_SW_DWRkw_T_-15_-10.png',
              inverty=inverty, figax=None,
              bins=(x,y),density=True, CFAD=True)

plt.close('all')