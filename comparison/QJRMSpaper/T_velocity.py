#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 19:48:03 2019

@author: dori
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('..')
from READ import slice_data
from READ import read_variables
from statistic import hist_and_plot

accMins = 5
freq = str(accMins)+'min'
iconfile = '../data/precipitation_icon.h5'
icon = pd.read_hdf(iconfile, key='stat')
icon = icon.reset_index().drop_duplicates(subset='index',
                                          keep='last').set_index('index')
pluviofile = '../data/precipitation_pluvio.h5'
pluvio = pd.read_hdf(pluviofile, key='stat')
mrrfile = '../data/precipitation_mrr.h5'
mrr = pd.read_hdf(mrrfile, key='stat')
mrr = mrr.resample(freq).apply(np.nansum)

ixi = pd.date_range(start='2015-11-11', end='2016-1-4', freq='9s')
icon.reindex(ixi)
icon = icon.resample(freq).apply(np.nansum)

ixp = pd.date_range(start='2015-11-11', end='2016-1-4', freq='1min')
pluvio.reindex(ixp)
pluvio = pluvio.resample(freq).apply(np.nansum)

pamtra = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='all_hydro', suffix='pamtra_icon.h5', pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T', 'unixtime',
                                 'V10', 'V35', 'V94',
                                 'W10', 'W35', 'W94'], minhour=6.0)

ice = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='only_ice', suffix='pamtra_icon.h5', pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T',
                                 'V10', 'V35', 'V94', 'unixtime',
                                 'W10', 'W35', 'W94'], minhour=6.0)

snow = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='only_snow', suffix='pamtra_icon.h5', pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T',
                                 'V10', 'V35', 'V94', 'unixtime',
                                 'W10', 'W35', 'W94'], minhour=6.0)

radar = read_variables(path='/work/develop/pamtraICON/comparison/data/radar/',
                       hydroset='', suffix='radar_regrid.h5', minhour=6.0,
                       varlist=['Z10', 'Z35', 'Z94', 'T',
                                'V10m5', 'V35m5', 'V94m5',
                                'W10', 'W35', 'W94', 'unixtime',
                                'quality_x', 'quality_w'])

radar.unixtime = pd.to_datetime(radar.unixtime.astype(np.int64), unit='s')
pamtra.unixtime = pd.to_datetime(pamtra.unixtime.astype(np.int64), unit='s')

radar['RR'] = (pluvio.resample('1s').nearest().loc[radar.unixtime]*60/accMins).values
pamtra['RR'] = (icon.resample('1s').nearest().loc[pamtra.unixtime]*60/accMins).values

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
                              yvar='T', xvar='V35m5',
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
                  yvar='T', xvar='V35m5',
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

## Low precipitation
minRR = -1.0
maxRR = 1.0
pre = 'LOW'
f, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2, figsize=(10.5, 9.))
r = hist_and_plot(slice_data(pamtra, 'RR', minvalue=minRR, maxvalue=maxRR),
                  'Simulated MDV Ka',
                  yvar='T', xvar='V35',
                  xlabel='MDV   [m/s]', ylabel='T   [deg C]',
                  vminmax=[0.1, 30],
                  xlim=[-5, 1], ylim=[-30, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VS.png',
                  inverty=inverty, figax=(f, ax11), stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radar, 'RR', minvalue=minRR, maxvalue=maxRR),
                  'Measured MDV Ka',
                  yvar='T', xvar='V35m5',
                  xlabel='MDV   [m/s]', ylabel='T   [deg C]',
                  vminmax=[0.1, 30],
                  xlim=[-5, 1], ylim=[-30, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VS.png',
                  inverty=inverty, figax=(f, ax12), stats=stats,
                  bins=(r[4], r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radar, 'RR', minvalue=minRR, maxvalue=maxRR),
                  'Measured SW Ka',
                  yvar='T', xvar='W35',
                  xlabel='SW   [m/s]',
                  ylabel='T   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 3], ylim=[-30, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VS.png',
                  inverty=inverty, figax=(f, ax22), stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(pamtra, 'RR', minvalue=minRR, maxvalue=maxRR),
                  'Simulated SW Ka',
                  yvar='T', xvar='W35',
                  xlabel='SW   [m/s]',
                  ylabel='T   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 3], ylim=[-30, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VS.png',
                  inverty=inverty, figax=(f, ax21), stats=stats,
                  bins=(r[4], r[5]), density=density, CFAD=CFAD)


f.suptitle('T-MDV  CFADs RR<1 mm/h', fontsize=12, fontweight='heavy', y=0.99)
f.tight_layout(pad=1.5, h_pad=0.5, w_pad=0.5)
f.text(x=0.5, y=0.49, s='T-SW   CFADs', fontsize=12, fontweight='heavy',
       horizontalalignment='center')
f.savefig(pre+'pamRad_T_VS.png', dpi=300)

## Medium precipitation
#minRR = 0.5
#maxRR = 1.0
#pre = 'MID'
#f, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2, figsize=(10.5, 9.))
#r = hist_and_plot(slice_data(pamtra, 'RR', minvalue=minRR, maxvalue=maxRR),
#                  'Simulated MDV Ka',
#                  yvar='T', xvar='V35',
#                  xlabel='MDV   [m/s]', ylabel='T   [deg C]',
#                  vminmax=[0.1, 30],
#                  xlim=[-5, 1], ylim=[-30, 10], lognorm=logrule,
#                  savename=pre+'pamRad_T_VS.png',
#                  inverty=inverty, figax=(f, ax11), stats=stats,
#                  bins=bins, density=density, CFAD=CFAD)
#
#r = hist_and_plot(slice_data(radar, 'RR', minvalue=minRR, maxvalue=maxRR),
#                  'Measured MDV Ka',
#                  yvar='T', xvar='V35m5',
#                  xlabel='MDV   [m/s]', ylabel='T   [deg C]',
#                  vminmax=[0.1, 30],
#                  xlim=[-5, 1], ylim=[-30, 10], lognorm=logrule,
#                  savename=pre+'pamRad_T_VS.png',
#                  inverty=inverty, figax=(f, ax12), stats=stats,
#                  bins=(r[4], r[5]), density=density, CFAD=CFAD)
#
#r = hist_and_plot(slice_data(radar, 'RR', minvalue=minRR, maxvalue=maxRR),
#                  'Measured SW Ka',
#                  yvar='T', xvar='W35',
#                  xlabel='SW   [m/s]',
#                  ylabel='T   [deg C]',
#                  vminmax=[0.1, 40],
#                  xlim=[0, 3], ylim=[-30, 10], lognorm=logrule,
#                  savename=pre+'pamRad_T_VS.png',
#                  inverty=inverty, figax=(f, ax22), stats=stats,
#                  bins=bins, density=density, CFAD=CFAD)
#
#r = hist_and_plot(slice_data(pamtra, 'RR', minvalue=minRR, maxvalue=maxRR),
#                  'Simulated SW Ka',
#                  yvar='T', xvar='W35',
#                  xlabel='SW   [m/s]',
#                  ylabel='T   [deg C]',
#                  vminmax=[0.1, 40],
#                  xlim=[0, 3], ylim=[-30, 10], lognorm=logrule,
#                  savename=pre+'pamRad_T_VS.png',
#                  inverty=inverty, figax=(f, ax21), stats=stats,
#                  bins=(r[4], r[5]), density=density, CFAD=CFAD)
#
#
#f.suptitle('T-MDV  CFADs 0.5<RR<1 mm/h', fontsize=12, fontweight='heavy', y=0.99)
#f.tight_layout(pad=1.5, h_pad=0.5, w_pad=0.5)
#f.text(x=0.5, y=0.49, s='T-SW   CFADs', fontsize=12, fontweight='heavy',
#       horizontalalignment='center')
#f.savefig(pre+'pamRad_T_VS.png', dpi=300)

## High precipitation
minRR = 1.0
maxRR = 91.0
pre = 'HIG'
f, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2, figsize=(10.5, 9.))
r = hist_and_plot(slice_data(pamtra, 'RR', minvalue=minRR, left=True),
                  'Simulated MDV Ka',
                  yvar='T', xvar='V35',
                  xlabel='MDV   [m/s]', ylabel='T   [deg C]',
                  vminmax=[0.1, 30],
                  xlim=[-5, 1], ylim=[-30, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VS.png',
                  inverty=inverty, figax=(f, ax11), stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radar, 'RR', minvalue=minRR, left=True),
                  'Measured MDV Ka',
                  yvar='T', xvar='V35m5',
                  xlabel='MDV   [m/s]', ylabel='T   [deg C]',
                  vminmax=[0.1, 30],
                  xlim=[-5, 1], ylim=[-30, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VS.png',
                  inverty=inverty, figax=(f, ax12), stats=stats,
                  bins=(r[4], r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radar, 'RR', minvalue=minRR, left=True),
                  'Measured SW Ka',
                  yvar='T', xvar='W35',
                  xlabel='SW   [m/s]',
                  ylabel='T   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 3], ylim=[-30, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VS.png',
                  inverty=inverty, figax=(f, ax22), stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(pamtra, 'RR', minvalue=minRR, left=True),
                  'Simulated SW Ka',
                  yvar='T', xvar='W35',
                  xlabel='SW   [m/s]',
                  ylabel='T   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 3], ylim=[-30, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VS.png',
                  inverty=inverty, figax=(f, ax21), stats=stats,
                  bins=(r[4], r[5]), density=density, CFAD=CFAD)


f.suptitle('T-MDV  CFADs RR>1 mm/h', fontsize=12, fontweight='heavy', y=0.99)
f.tight_layout(pad=1.5, h_pad=0.5, w_pad=0.5)
f.text(x=0.5, y=0.49, s='T-SW   CFADs', fontsize=12, fontweight='heavy',
       horizontalalignment='center')
f.savefig(pre+'pamRad_T_VS.png', dpi=300)