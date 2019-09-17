#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 11:49:02 2019

@author: dori
"""

import sys
sys.path.append('..')
from READ import slice_data
from READ import read_variables
from statistic import hist_and_plot
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import time as timer

pamtra = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='all_hydro', suffix='pamtra_icon.h5',
                        minhour=6.0, pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T', 'Hgt', 'unixtime'])


radar = read_variables(path='/work/develop/pamtraICON/comparison/data/radar/',
                       hydroset='', suffix='radar_regrid.h5', minhour=6.0,
                       varlist=['Z10', 'Z35', 'Z94', 'T', 'Hgt',
                                'quality_x', 'quality_w', 'unixtime']
                      ).reset_index()

radpure = read_variables(path='/work/develop/pamtraICON/comparison/data/radar/',
                       hydroset='', suffix='radar.h5', minhour=6.0,
                       varlist=['Z10', 'Z35', 'Z94', 'T', 'Hgt',
                                'quality_x', 'quality_w', 'unixtime']
                      ).reset_index()

#radarw = slice_data(radar, 'quality_w', maxvalue=8192)
#radarx = slice_data(radar, 'quality_x', maxvalue=8192)

pamtraProfiles = pamtra.dropna(subset=['Z35']).groupby('unixtime')
#radarxProfiles = radarx.groupby('unixtime')
#radarwProfiles = radarw.groupby('unixtime')
radarProfiles = radar.dropna(subset=['Z35']).groupby('groupT')
radpuProfiles = radpure.dropna(subset=['Z35']).groupby('unixtime')

def reduction(x):
  d = {}
  y = x.loc[x['T'] < 0]
  try:
    CTidx = y[['Hgt','Z35']].dropna().sort_values('Hgt').index[-1]
    d['CTT'] = y.loc[CTidx]['T']
    d['CTH'] = y.loc[CTidx]['Hgt']
  except IndexError:
    d['CTT'] = np.nan
    d['CTH'] = np.nan
  try:
    yw = y.loc[y['quality_w'] < 8192]
  except KeyError:
    yw = y
  try:
    kwIdx = yw[['DWRkw']].dropna().sort_values('DWRkw').index[-1]
    d['maxDWRkw'] = yw.loc[kwIdx]['DWRkw']
    d['maxDWRkwH'] = yw.loc[kwIdx]['Hgt']
  except IndexError:
    d['maxDWRkw'] = np.nan
    d['maxDWRkwH'] = np.nan
  try:
    yx = y.loc[y['quality_x'] < 8192]
  except KeyError:
    yx = y
  try:
    xkIdx = yx[['DWRxk']].dropna().sort_values('DWRxk').index[-1]
    d['maxDWRxk'] = yx.loc[xkIdx]['DWRkw']
    d['maxDWRxkH'] = yx.loc[xkIdx]['Hgt']
  except IndexError:
    d['maxDWRxk'] = np.nan
    d['maxDWRxkH'] = np.nan
  return pd.Series(d, index=d.keys())

start = timer.time()
radpuStats = radpuProfiles.apply(reduction)
end = timer.time()
print(end-start)

start = timer.time()
radarStats = radarProfiles.apply(reduction)
end = timer.time()
print(end - start)
pamtraStats = pamtraProfiles.apply(reduction)
print(timer.time()-end)

rad = radarStats.dropna(how='all')
pam = pamtraStats.dropna(how='all')
rpu = radpuStats.dropna(how='all')

rad.to_hdf('../data/radStatCTTmaxDWRregrid.h5', key='stat')
pam.to_hdf('../data/pamStatCTTmaxDWR.h5', key='stat')
rpu.to_hdf('../data/radStatCTTmaxDWR.h5', key='stat')

logrule = True
density = False
CFAD = True
inverty = True
bins = 100
stats = ['mean', 'median', 'quantile']

r = hist_and_plot(pam, 'Simulated Cloud Top Temperature - max DWRkw',
                  yvar='CTT', xvar='maxDWRkw',
                  xlabel='max DWRkw   [dB]', ylabel='CTT   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 20], ylim=[-60, -10], lognorm=logrule,
                  savename='pam_CTT_DWRkw.png',
                  inverty=inverty, figax=None, stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(rad, 'Measured Regridded Cloud Top Temperature - max DWRkw',
                  yvar='CTT', xvar='maxDWRkw',
                  xlabel='max DWRkw   [dB]', ylabel='CTT   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 20], ylim=[-60, -10], lognorm=logrule,
                  savename='rad_CTT_DWRkw.png',
                  inverty=inverty, figax=None, stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(rpu, 'Measured Cloud Top Temperature - max DWRkw',
                  yvar='CTT', xvar='maxDWRkw',
                  xlabel='max DWRkw   [dB]', ylabel='CTT   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 20], ylim=[-60, -10], lognorm=logrule,
                  savename='rpu_CTT_DWRkw.png',
                  inverty=inverty, figax=None, stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(pam, 'Simulated Cloud Top Temperature - max DWRxk',
                  yvar='CTT', xvar='maxDWRxk',
                  xlabel='max DWRxk   [dB]', ylabel='CTT   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 20], ylim=[-60, -10], lognorm=logrule,
                  savename='pam_CTT_DWRxk.png',
                  inverty=inverty, figax=None, stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(rad, 'Measured Regridded Cloud Top Temperature - max DWRxk',
                  yvar='CTT', xvar='maxDWRxk',
                  xlabel='max DWRxk   [dB]', ylabel='CTT   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 20], ylim=[-60, -10], lognorm=logrule,
                  savename='rad_CTT_DWRxk.png',
                  inverty=inverty, figax=None, stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(rpu, 'Measured Cloud Top Temperature - max DWRxk',
                  yvar='CTT', xvar='maxDWRxk',
                  xlabel='max DWRxk   [dB]', ylabel='CTT   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 20], ylim=[-60, -10], lognorm=logrule,
                  savename='rpu_CTT_DWRxk.png',
                  inverty=inverty, figax=None, stats=stats,
                  bins=bins, density=density, CFAD=CFAD)