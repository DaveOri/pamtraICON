#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 19:34:35 2019

@author: dori
"""

import sys
sys.path.append('..')
from READ import slice_data
from READ import read_variables
from statistic import hist_and_plot
import matplotlib.pyplot as plt
import netCDF4
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

model = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='all_hydro', suffix='pamtra_icon.h5', pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T',
                                 'V10', 'V35', 'V94'], minhour=6.0)

#ice = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
#                        hydroset='only_ice', suffix='pamtra_icon.h5', pamtra=True,
#                        varlist=['Z10', 'Z35', 'Z94', 'T',
#                                 'V10', 'V35', 'V94'], minhour=6.0)
#
#snow = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
#                        hydroset='only_snow', suffix='pamtra_icon.h5', pamtra=True,
#                        varlist=['Z10', 'Z35', 'Z94', 'T',
#                                 'V10', 'V35', 'V94'], minhour=6.0)

radar = read_variables(path='/work/develop/pamtraICON/comparison/data/radar/',
                       hydroset='', suffix='radar.h5', minhour=6.0,
                       varlist=['Z10', 'Z35', 'Z94', 'T',
                                'V10avg', 'V35avg', 'V94avg', 'quality_x', 'quality_w'])

radarw = slice_data(radar, 'quality_w', maxvalue=8192)
radarx = slice_data(radar, 'quality_x', maxvalue=8192)

data = netCDF4.Dataset('../data/idealized_hydro.nc')
data_simple = netCDF4.Dataset('../data/idealized_hydro_simple.nc')
datavar = data.variables
Zvar = data_simple.variables['Ze']
Ze = Zvar[:]
Ze[:,0,:,2,0,0] = Ze[:,0,:,2,0,0] #+ 10.*np.log10(0.93/0.72)
V = datavar['Radar_MeanDopplerVel'][:]


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

Zc[:,0] -= 0.5
Zr[:,0] -= 0.5

logrule = True
density = True
CFAD = False
inverty = True
bins = 100
lw = 3
vminmax = [1e2, 1e6]
pamtra = slice_data(model, 'Z10', -14)
r = hist_and_plot(slice_data(pamtra, 'T', -20, -5),
                              'Simulated -20<T<-5',
                              yvar='V35', xvar='DWRxk',
                              xlabel='DWRxk   [dB]', ylabel='MDV   [m/s]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-3, 0], lognorm=logrule,
                              savename='pamtra_Vk_DWRxk.png',
                              inverty=inverty, figax=None,
                              bins=bins, density=density, CFAD=CFAD)
plt.gca().set_prop_cycle(color=colors[1:])
plt.plot(Zi[:,0]-Zi[:,1], -Vi[:,0], label='ice crystals', lw=lw)
plt.plot(Zr[:,0]-Zr[:,1], -Vr[:,0], label='raindrops', lw=lw)
plt.plot(Zs[:,0]-Zs[:,1], -Vs[:,0], label='snowflakes', lw=lw)
plt.plot(Zg[:,0]-Zg[:,1], -Vg[:,0], label='graupel', lw=lw)
plt.plot(Zh[:,0]-Zh[:,1], -Vh[:,0], label='hail', lw=lw)
plt.savefig('pamtra_Vk_DWRxk.png', dpi=300)

r = hist_and_plot(slice_data(radarx, 'T', -20, -5),
                              'Measured -20<T<-5',
                              yvar='V35avg', xvar='DWRxk',
                              xlabel='DWRxk   [dB]', ylabel='MDV   [m/s]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-3, 0], lognorm=logrule,
                              savename='radar_Vk_DWRxk.png',
                              inverty=inverty, figax=None,
                              bins=bins, density=density, CFAD=CFAD)
plt.gca().set_prop_cycle(color=colors[1:])
plt.plot(Zi[:,0]-Zi[:,1], -Vi[:,0], label='ice crystals', lw=lw)
plt.plot(Zr[:,0]-Zr[:,1], -Vr[:,0], label='raindrops', lw=lw)
plt.plot(Zs[:,0]-Zs[:,1], -Vs[:,0], label='snowflakes', lw=lw)
plt.plot(Zg[:,0]-Zg[:,1], -Vg[:,0], label='graupel', lw=lw)
plt.plot(Zh[:,0]-Zh[:,1], -Vh[:,0], label='hail', lw=lw)
plt.legend(loc=1)
plt.savefig('radar_Vk_DWRxk.png', dpi=300)

r = hist_and_plot(slice_data(pamtra, 'T', -20, -5),
                              'Simulated -20<T<-5',
                              yvar='V35', xvar='DWRkw',
                              xlabel='DWRkw   [dB]', ylabel='MDV   [m/s]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-3, 0], lognorm=logrule,
                              savename='pamtra_Vk_DWRkw.png',
                              inverty=inverty, figax=None,
                              bins=bins, density=density, CFAD=CFAD)
plt.gca().set_prop_cycle(color=colors[1:])
plt.plot(Zi[:,1]-Zi[:,2], -Vi[:,0], label='ice crystals', lw=lw)
plt.plot(Zr[:,1]-Zr[:,2], -Vr[:,0], label='raindrops', lw=lw)
plt.plot(Zs[:,1]-Zs[:,2], -Vs[:,0], label='snowflakes', lw=lw)
plt.plot(Zg[:,1]-Zg[:,2], -Vg[:,0], label='graupel', lw=lw)
plt.plot(Zh[:,1]-Zh[:,2], -Vh[:,0], label='hail', lw=lw)
plt.savefig('pamtra_Vk_DWRkw.png', dpi=300)

r = hist_and_plot(slice_data(radarw, 'T', -20, -5),
                              'Measured -20<T<-5',
                              yvar='V35avg', xvar='DWRkw',
                              xlabel='DWRkw   [dB]', ylabel='MDV   [m/s]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-3, 0], lognorm=logrule,
                              savename='radar_Vk_DWRkw.png',
                              inverty=inverty, figax=None,
                              bins=bins, density=density, CFAD=CFAD)
plt.gca().set_prop_cycle(color=colors[1:])
plt.plot(Zi[:,1]-Zi[:,2], -Vi[:,0], label='ice crystals', lw=lw)
plt.plot(Zr[:,1]-Zr[:,2], -Vr[:,0], label='raindrops', lw=lw)
plt.plot(Zs[:,1]-Zs[:,2], -Vs[:,0], label='snowflakes', lw=lw)
plt.plot(Zg[:,1]-Zg[:,2], -Vg[:,0], label='graupel', lw=lw)
plt.plot(Zh[:,1]-Zh[:,2], -Vh[:,0], label='hail', lw=lw)
plt.legend(loc=1)
plt.savefig('radar_Vk_DWRkw.png', dpi=300)

#%%

r = hist_and_plot(slice_data(pamtra, 'T', -5, 40),
                              'Simulated -5<T',
                              yvar='V35', xvar='DWRxk',
                              xlabel='DWRxk   [dB]', ylabel='MDV   [m/s]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-3, 0], lognorm=logrule,
                              savename='pamtra_Vk_DWRxk_T>-5.png',
                              inverty=inverty, figax=None,
                              bins=bins, density=density, CFAD=CFAD)
plt.gca().set_prop_cycle(color=colors[1:])
plt.plot(Zi[:,0]-Zi[:,1], -Vi[:,0], label='ice crystals', lw=lw)
plt.plot(Zr[:,0]-Zr[:,1], -Vr[:,0], label='raindrops', lw=lw)
plt.plot(Zs[:,0]-Zs[:,1], -Vs[:,0], label='snowflakes', lw=lw)
plt.plot(Zg[:,0]-Zg[:,1], -Vg[:,0], label='graupel', lw=lw)
plt.plot(Zh[:,0]-Zh[:,1], -Vh[:,0], label='hail', lw=lw)
plt.savefig('pamtra_Vk_DWRxk_T>-5.png', dpi=300)

r = hist_and_plot(slice_data(radarx, 'T', -5, 40),
                              'Measured -5<T',
                              yvar='V35avg', xvar='DWRxk',
                              xlabel='DWRxk   [dB]', ylabel='MDV   [m/s]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-3, 0], lognorm=logrule,
                              savename='radar_Vk_DWRxk_T>-5.png',
                              inverty=inverty, figax=None,
                              bins=bins, density=density, CFAD=CFAD)
plt.gca().set_prop_cycle(color=colors[1:])
plt.plot(Zi[:,0]-Zi[:,1], -Vi[:,0], label='ice crystals', lw=lw)
plt.plot(Zr[:,0]-Zr[:,1], -Vr[:,0], label='raindrops', lw=lw)
plt.plot(Zs[:,0]-Zs[:,1], -Vs[:,0], label='snowflakes', lw=lw)
plt.plot(Zg[:,0]-Zg[:,1], -Vg[:,0], label='graupel', lw=lw)
plt.plot(Zh[:,0]-Zh[:,1], -Vh[:,0], label='hail', lw=lw)
plt.legend(loc=1)
plt.savefig('radar_Vk_DWRxk_T>-5.png', dpi=300)

r = hist_and_plot(slice_data(pamtra, 'T', -5, 40),
                              'Simulated -5<T',
                              yvar='V35', xvar='DWRkw',
                              xlabel='DWRkw   [dB]', ylabel='MDV   [m/s]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-3, 0], lognorm=logrule,
                              savename='pamtra_Vk_DWRkw_T>-5.png',
                              inverty=inverty, figax=None,
                              bins=bins, density=density, CFAD=CFAD)
plt.gca().set_prop_cycle(color=colors[1:])
plt.plot(Zi[:,1]-Zi[:,2], -Vi[:,0], label='ice crystals', lw=lw)
plt.plot(Zr[:,1]-Zr[:,2], -Vr[:,0], label='raindrops', lw=lw)
plt.plot(Zs[:,1]-Zs[:,2], -Vs[:,0], label='snowflakes', lw=lw)
plt.plot(Zg[:,1]-Zg[:,2], -Vg[:,0], label='graupel', lw=lw)
plt.plot(Zh[:,1]-Zh[:,2], -Vh[:,0], label='hail', lw=lw)
plt.savefig('pamtra_Vk_DWRkw_T>-5.png', dpi=300)

r = hist_and_plot(slice_data(radarw, 'T', -5, 40),
                              'Measured -5<T',
                              yvar='V35avg', xvar='DWRkw',
                              xlabel='DWRkw   [dB]', ylabel='MDV   [m/s]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-3, 0], lognorm=logrule,
                              savename='radar_Vk_DWRkw_T>-5.png',
                              inverty=inverty, figax=None,
                              bins=bins, density=density, CFAD=CFAD)
plt.gca().set_prop_cycle(color=colors[1:])
plt.plot(Zi[:,1]-Zi[:,2], -Vi[:,0], label='ice crystals', lw=lw)
plt.plot(Zr[:,1]-Zr[:,2], -Vr[:,0], label='raindrops', lw=lw)
plt.plot(Zs[:,1]-Zs[:,2], -Vs[:,0], label='snowflakes', lw=lw)
plt.plot(Zg[:,1]-Zg[:,2], -Vg[:,0], label='graupel', lw=lw)
plt.plot(Zh[:,1]-Zh[:,2], -Vh[:,0], label='hail', lw=lw)
plt.legend(loc=1)
plt.savefig('radar_Vk_DWRkw_T>-5.png', dpi=300)
