#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 18:53:47 2019

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

rad = pd.read_hdf('../data/radStatCTTmaxDWRregrid.h5', key='stat')
rad.index = pd.to_datetime(rad.index, unit='s')
rad = rad.resample('1s').nearest(limit=1)
rpu = pd.read_hdf('../data/radStatCTTmaxDWR.h5', key='stat')
rpu.index = pd.to_datetime(rpu.index, unit='s')
rpu = rpu.resample('1s').nearest(limit=1)
pam = pd.read_hdf('../data/pamStatCTTmaxDWR.h5', key='stat')
pam.index = pd.to_datetime(pam.index, unit='s')
pam = pam.resample('1s').nearest(limit=1)

ixi = pd.date_range(start='2015-11-11', end='2016-1-4', freq='9s')
icon.reindex(ixi)
icon = icon.resample(freq).apply(np.nansum)

ixp = pd.date_range(start='2015-11-11', end='2016-1-4', freq='1min')
pluvio.reindex(ixp)
pluvio = pluvio.resample(freq).apply(np.nansum)

pamtra = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='all_hydro', suffix='pamtra_icon.h5',
                        pamtra=True, minhour=6.0, attHYD=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T', 'unixtime',
                                 'H10', 'H35', 'H94'])

radar = read_variables(path='/work/develop/pamtraICON/comparison/data/radar/',
                       hydroset='', suffix='radar_regrid.h5', minhour=6.0,
                       varlist=['Z10', 'Z35', 'Z94', 'T', 'unixtime',
                                'quality_x', 'quality_w'])

radar = read_variables(path='/work/develop/pamtraICON/comparison/data/radar/',
                       hydroset='', suffix='radar.h5', minhour=6.0,
                       varlist=['Z10', 'Z35', 'Z94', 'T', 'unixtime',
                                'quality_x', 'quality_w'])

radar.unixtime = pd.to_datetime(radar.unixtime.astype(np.int64), unit='s')
pamtra.unixtime = pd.to_datetime(pamtra.unixtime.astype(np.int64), unit='s')

radar['RR'] = (pluvio.resample('1s').nearest().loc[radar.unixtime]*60/accMins).values
pamtra['RR'] = (icon.resample('1s').nearest().loc[pamtra.unixtime]*60/accMins).values
rad = rad.loc[radar.unixtime]
rpu = rpu.loc[radar.unixtime]
pam = pam.loc[pamtra.unixtime]

#for col in rpu.columns:
#  radar[col] = rpu[col].values
for col in rad.columns:
  radar[col] = rad[col].values
for col in pam.columns:
  pamtra[col] = pam[col].values

pamtra = slice_data(pamtra, 'Z10', minvalue=-7)
#radarw = slice_data(radar, 'quality_w', maxvalue=8192)
#radarx = slice_data(radar, 'quality_x', maxvalue=8192)
radarw = slice_data(radar, 'quality_w', maxvalue=4096)
radarx = slice_data(radar, 'quality_x', maxvalue=4096)

logrule = True
density = True
CFAD = True
bins = (np.arange(-5, 15, 0.35), np.arange(-45,1,0.6))
vminmax = [0.1, 20]
stats = ['mean', 'median', 'quartile', 'decile']

r = hist_and_plot(pamtra, 'Simulated CFAD T - DWRxk',
                              yvar='T', xvar='DWRxk',
                              xlabel='DWRxk   [dB]', ylabel='T   [deg C]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                              savename='pamtra_T_DWRxk.png',
                              inverty=True, figax=None, stats=stats,
                              bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(radarx, 'Measured CFAD T - DWRxk',
                              yvar='T', xvar='DWRxk',
                              xlabel='DWRxk   [dB]', ylabel='T   [deg C]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                              savename='radar_T_DWRxk.png',
                              inverty=True, figax=None, stats=stats,
                              bins=(r[4],r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(pamtra, 'Simulated CFAD T - DWRkw',
                              yvar='T', xvar='DWRkw',
                              xlabel='DWRkw   [dB]', ylabel='T   [deg C]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                              savename='pamtra_T_DWRkw.png',
                              inverty=True, figax=None, stats=stats,
                              bins=(r[4],r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(radarw, 'Measured CFAD T - DWRkw',
                              yvar='T', xvar='DWRkw',
                              xlabel='DWRkw   [dB]', ylabel='T   [deg C]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                              savename='radar_T_DWRkw.png',
                              inverty=True, figax=None, stats=stats,
                              bins=(r[4],r[5]), density=density, CFAD=CFAD)

#%% Combined plot for paper

f, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2, figsize=(10.5, 9.))
r = hist_and_plot(pamtra, 'Simulated Ka-W',
                              yvar='T', xvar='DWRkw',
                              xlabel='DWR$_{K_aW}$   [dB]', ylabel='T   [deg C]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                              savename='pamRad_T_DWRs.png',
                              inverty=True, figax=(f, ax11), stats=stats,
                              bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(radarw, 'Measured Ka-W',
                              yvar='T', xvar='DWRkw',
                              xlabel='DWR$_{K_aW}$   [dB]', ylabel='T   [deg C]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                              savename='pamRad_T_DWRs.png',
                              inverty=True, figax=(f, ax12), stats=stats,
                              bins=(r[4],r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(pamtra, 'Simulated X-Ka',
                              yvar='T', xvar='DWRxk',
                              xlabel='DWR$_{XK_a}$   [dB]', ylabel='T   [deg C]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                              savename='pamRad_T_DWRs.png',
                              inverty=True, figax=(f, ax21), stats=stats,
                              bins=(r[4],r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(radarx, 'Measured X-Ka',
                              yvar='T', xvar='DWRxk',
                              xlabel='DWR$_{XK_a}$   [dB]', ylabel='T   [deg C]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                              savename='pamRad_T_DWRs.png',
                              inverty=True, figax=(f, ax22), stats=stats,
                              bins=(r[4],r[5]), density=density, CFAD=CFAD)

f.suptitle('T-DWRs CFADs', fontsize=12, fontweight='heavy', y=0.99)
f.tight_layout(pad=1.5, h_pad=0.5, w_pad=0.5)
f.savefig('pamRad_T_DWRs.pdf', dpi=300)
f.savefig('pamRad_T_DWRs.png', dpi=300)

#%% Low precipitation
maxRR = 1.0
pre = 'LOW'
f, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2, figsize=(10.5, 9.))
r = hist_and_plot(slice_data(pamtra, 'RR', maxvalue=maxRR),
                  'Simulated Ka-W',
                  yvar='T', xvar='DWRkw',
                  xlabel='DWR$_{K_aW}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax11), stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radarw, 'RR', maxvalue=maxRR),
                  'Measured Ka-W',
                  yvar='T', xvar='DWRkw',
                  xlabel='DWR$_{K_aW}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax12), stats=stats,
                  bins=(r[4],r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(pamtra, 'RR', maxvalue=maxRR),
                  'Simulated X-Ka',
                  yvar='T', xvar='DWRxk',
                  xlabel='DWR$_{XK_a}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax21), stats=stats,
                  bins=(r[4],r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radarx, 'RR', maxvalue=maxRR),
                  'Measured X-Ka',
                  yvar='T', xvar='DWRxk',
                  xlabel='DWR$_{XK_a}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax22), stats=stats,
                  bins=(r[4],r[5]), density=density, CFAD=CFAD)

f.suptitle('T-DWRs CFADs RR<1 mm/h', fontsize=12, fontweight='heavy', y=0.99)
f.tight_layout(pad=1.5, h_pad=0.5, w_pad=0.5)
f.savefig(pre+'pamRad_T_DWRs.pdf', dpi=300)
f.savefig(pre+'pamRad_T_DWRs.png', dpi=300)

##%% Mid precipitation
#minRR = 0.5
#maxRR = 1.0
#pre = 'MID'
#f, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2, figsize=(10.5, 9.))
#r = hist_and_plot(slice_data(pamtra, 'RR', maxvalue=maxRR),
#                  'Simulated Ka-W',
#                  yvar='T', xvar='DWRkw',
#                  xlabel='DWR$_{K_aW}$   [dB]', ylabel='T   [deg C]',
#                  vminmax=vminmax,
#                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
#                  savename=pre+'pamRad_T_DWRs.png',
#                  inverty=True, figax=(f, ax11), stats=stats,
#                  bins=bins, density=density, CFAD=CFAD)
#
#r = hist_and_plot(slice_data(radarw, 'RR', maxvalue=maxRR),
#                  'Measured Ka-W',
#                  yvar='T', xvar='DWRkw',
#                  xlabel='DWR$_{K_aW}$   [dB]', ylabel='T   [deg C]',
#                  vminmax=vminmax,
#                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
#                  savename=pre+'pamRad_T_DWRs.png',
#                  inverty=True, figax=(f, ax12), stats=stats,
#                  bins=(r[4],r[5]), density=density, CFAD=CFAD)
#
#r = hist_and_plot(slice_data(pamtra, 'RR', maxvalue=maxRR),
#                  'Simulated X-Ka',
#                  yvar='T', xvar='DWRxk',
#                  xlabel='DWR$_{XK_a}$   [dB]', ylabel='T   [deg C]',
#                  vminmax=vminmax,
#                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
#                  savename=pre+'pamRad_T_DWRs.png',
#                  inverty=True, figax=(f, ax21), stats=stats,
#                  bins=(r[4],r[5]), density=density, CFAD=CFAD)
#
#r = hist_and_plot(slice_data(radarx, 'RR', maxvalue=maxRR),
#                  'Measured X-Ka',
#                  yvar='T', xvar='DWRxk',
#                  xlabel='DWR$_{XK_a}$   [dB]', ylabel='T   [deg C]',
#                  vminmax=vminmax,
#                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
#                  savename=pre+'pamRad_T_DWRs.png',
#                  inverty=True, figax=(f, ax22), stats=stats,
#                  bins=(r[4],r[5]), density=density, CFAD=CFAD)
#
#f.suptitle('T-DWRs CFADs 0.5<RR<1 mm/h', fontsize=12, fontweight='heavy', y=0.99)
#f.tight_layout(pad=1.5, h_pad=0.5, w_pad=0.5)
#f.savefig(pre+'pamRad_T_DWRs.png', dpi=300)

#%% High precipitation
minRR = 1.0
pre = 'HIG'
f, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2, figsize=(10.5, 9.))
r = hist_and_plot(slice_data(pamtra, 'RR', minvalue=minRR, left=True),
                  'Simulated Ka-W',
                  yvar='T', xvar='DWRkw',
                  xlabel='DWR$_{K_aW}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax11), stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radarw, 'RR', minvalue=minRR, left=True),
                  'Measured Ka-W',
                  yvar='T', xvar='DWRkw',
                  xlabel='DWR$_{K_aW}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax12), stats=stats,
                  bins=(r[4],r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(pamtra, 'RR', minvalue=minRR, left=True),
                  'Simulated X-Ka',
                  yvar='T', xvar='DWRxk',
                  xlabel='DWR$_{XK_a}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax21), stats=stats,
                  bins=(r[4],r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radarx, 'RR', minvalue=minRR, left=True),
                  'Measured X-Ka',
                  yvar='T', xvar='DWRxk',
                  xlabel='DWR$_{XK_a}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax22), stats=stats,
                  bins=(r[4],r[5]), density=density, CFAD=CFAD)

f.suptitle('T-DWRs CFADs RR>1 mm/h', fontsize=12, fontweight='heavy', y=0.99)
f.tight_layout(pad=1.5, h_pad=0.5, w_pad=0.5)
f.savefig(pre+'pamRad_T_DWRs.pdf', dpi=300)
f.savefig(pre+'pamRad_T_DWRs.png', dpi=300)

#%% High DWRkw
minKW = 6.0
pre = 'higKW'
f, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2, figsize=(10.5, 9.))
r = hist_and_plot(slice_data(pamtra, 'maxDWRkw', minvalue=minKW, left=True),
                  'Simulated Ka-W',
                  yvar='T', xvar='DWRkw',
                  xlabel='DWR$_{K_aW}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax11), stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radarw, 'maxDWRkw', minvalue=minKW, left=True),
                  'Measured Ka-W',
                  yvar='T', xvar='DWRkw',
                  xlabel='DWR$_{K_aW}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax12), stats=stats,
                  bins=(r[4],r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(pamtra, 'maxDWRkw', minvalue=minKW, left=True),
                  'Simulated X-Ka',
                  yvar='T', xvar='DWRxk',
                  xlabel='DWR$_{XK_a}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax21), stats=stats,
                  bins=(r[4],r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radarx, 'maxDWRkw', minvalue=minKW, left=True),
                  'Measured X-Ka',
                  yvar='T', xvar='DWRxk',
                  xlabel='DWR$_{XK_a}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax22), stats=stats,
                  bins=(r[4],r[5]), density=density, CFAD=CFAD)

f.suptitle('T-DWRs CFADs maxDWRkw>6 dB', fontsize=12, fontweight='heavy', y=0.99)
f.tight_layout(pad=1.5, h_pad=0.5, w_pad=0.5)
f.savefig(pre+'pamRad_T_DWRs.pdf', dpi=300)
f.savefig(pre+'pamRad_T_DWRs.png', dpi=300)

#%% Low DWRkw
minKW = 0.0
maxKW = 6.0
pre = 'lowKW'
f, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2, figsize=(10.5, 9.))
r = hist_and_plot(slice_data(pamtra, 'maxDWRkw',
                             minvalue=minKW, maxvalue=maxKW, left=True),
                  'Simulated Ka-W',
                  yvar='T', xvar='DWRkw',
                  xlabel='DWR$_{K_aW}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax11), stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radarw, 'maxDWRkw',
                             minvalue=minKW, maxvalue=maxKW, left=True),
                  'Measured Ka-W',
                  yvar='T', xvar='DWRkw',
                  xlabel='DWR$_{K_aW}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax12), stats=stats,
                  bins=(r[4],r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(pamtra, 'maxDWRkw',
                             minvalue=minKW, maxvalue=maxKW, left=True),
                  'Simulated X-Ka',
                  yvar='T', xvar='DWRxk',
                  xlabel='DWR$_{XK_a}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax21), stats=stats,
                  bins=(r[4],r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radarx, 'maxDWRkw',
                             minvalue=minKW, maxvalue=maxKW, left=True),
                  'Measured X-Ka',
                  yvar='T', xvar='DWRxk',
                  xlabel='DWR$_{XK_a}$   [dB]', ylabel='T   [deg C]',
                  vminmax=vminmax,
                  xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                  savename=pre+'pamRad_T_DWRs.png',
                  inverty=True, figax=(f, ax22), stats=stats,
                  bins=(r[4],r[5]), density=density, CFAD=CFAD)

f.suptitle('T-DWRs CFADs maxDWRkw<6 dB', fontsize=12, fontweight='heavy', y=0.99)
f.tight_layout(pad=1.5, h_pad=0.5, w_pad=0.5)
f.savefig(pre+'pamRad_T_DWRs.pdf', dpi=300)
f.savefig(pre+'pamRad_T_DWRs.png', dpi=300)