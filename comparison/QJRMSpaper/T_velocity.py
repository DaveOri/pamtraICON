#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 19:48:03 2019

@author: dori
"""

import sys
sys.path.append('..')
from READ import slice_data
from READ import read_variables
from statistic import hist_and_plot
import matplotlib.pyplot as plt
import netCDF4

pamtra = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='all_hydro', suffix='pamtra_icon.h5', pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T',
                                 'V10', 'V35', 'V94',
                                 'W10', 'W35', 'W94'], minhour=6.0)

ice = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='only_ice', suffix='pamtra_icon.h5', pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T',
                                 'V10', 'V35', 'V94',
                                 'W10', 'W35', 'W94'], minhour=6.0)

snow = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='only_snow', suffix='pamtra_icon.h5', pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T',
                                 'V10', 'V35', 'V94',
                                 'W10', 'W35', 'W94'], minhour=6.0)

radar = read_variables(path='/work/develop/pamtraICON/comparison/data/radar/',
                       hydroset='', suffix='radar_regrid.h5', minhour=6.0,
                       varlist=['Z10', 'Z35', 'Z94', 'T',
                                'V10avg', 'V35avg', 'V94avg',
                                'W10', 'W35', 'W94',
                                'quality_x', 'quality_w'])

radarw = slice_data(radar, 'quality_w', maxvalue=8192)
radarx = slice_data(radar, 'quality_x', maxvalue=8192)

logrule = True
density = False
CFAD = True
inverty = True
bins = 100
stats = ['mean', 'median', 'quartile', 'decile']

r = hist_and_plot(pamtra, 'Simulated CFAD   T - MDV',
                              yvar='T', xvar='V35',
                              xlabel='MDV ka   [m/s]', ylabel='T   [K]',
                              vminmax=[0.1, 30],
                              xlim=[-5, 1], ylim=[-30, 10], lognorm=logrule,
                              savename='pamtra_T_Vk.png',
                              inverty=inverty, figax=None, stats=stats,
                              bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(radar, 'Measured CFAD   T- MDV',
                              yvar='T', xvar='V35avg',
                              xlabel='MDV ka   [m/s]', ylabel='T   [K]',
                              vminmax=[0.1, 30],
                              xlim=[-5, 1], ylim=[-30, 10], lognorm=logrule,
                              savename='radar_T_Vk.png',
                              inverty=inverty, figax=None, stats=stats,
                              bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(pamtra, 'Simulated CFAD   T - SW',
                              yvar='T', xvar='W35',
                              xlabel='Spectral Width Ka   [m/s]', ylabel='T   [K]',
                              vminmax=[0.1, 40],
                              xlim=[0, 3], ylim=[-30, 10], lognorm=logrule,
                              savename='pamtra_T_SWk.png',
                              inverty=inverty, figax=None, stats=stats,
                              bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(radar, 'Measured CFAD   T- SW',
                              yvar='T', xvar='W35',
                              xlabel='Spectral Width Ka   [m/s]', ylabel='T   [K]',
                              vminmax=[0.1, 40],
                              xlim=[0, 3], ylim=[-30, 10], lognorm=logrule,
                              savename='radar_T_SWk.png',
                              inverty=inverty, figax=None, stats=stats,
                              bins=bins, density=density, CFAD=CFAD)
#%% Combined plot for paper

f, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2, figsize=(10.5, 9.))
r = hist_and_plot(pamtra, 'Simulated MDV Ka',
                  yvar='T', xvar='V35',
                  xlabel='MDV   [m/s]', ylabel='T   [deg C]',
                  vminmax=[0.1, 30],
                  xlim=[-5, 1], ylim=[-30, 10], lognorm=logrule,
                  savename='pamRad_T_VS.png',
                  inverty=inverty, figax=(f, ax11), stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(radar, 'Measured MDV Ka',
                  yvar='T', xvar='V35avg',
                  xlabel='MDV   [m/s]', ylabel='T   [deg C]',
                  vminmax=[0.1, 30],
                  xlim=[-5, 1], ylim=[-30, 10], lognorm=logrule,
                  savename='pamRad_T_VS.png',
                  inverty=inverty, figax=(f, ax12), stats=stats,
                  bins=(r[4], r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(radar, 'Measured SW Ka',
                  yvar='T', xvar='W35',
                  xlabel='SW   [m/s]',
                  ylabel='T   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 3], ylim=[-30, 10], lognorm=logrule,
                  savename='pamRad_T_VS.png',
                  inverty=inverty, figax=(f, ax22), stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(pamtra, 'Simulated SW Ka',
                  yvar='T', xvar='W35',
                  xlabel='SW   [m/s]',
                  ylabel='T   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 3], ylim=[-30, 10], lognorm=logrule,
                  savename='pamRad_T_VS.png',
                  inverty=inverty, figax=(f, ax21), stats=stats,
                  bins=(r[4], r[5]), density=density, CFAD=CFAD)


f.suptitle('T-MDV  CFADs', fontsize=12, fontweight='heavy', y=0.99)
f.tight_layout(pad=1.5, h_pad=0.5, w_pad=0.5)
f.text(x=0.5, y=0.49, s='T-SW   CFADs', fontsize=12, fontweight='heavy',
       horizontalalignment='center')
f.savefig('pamRad_T_VS.png', dpi=300)