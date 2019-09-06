#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 16:58:01 2019

@author: dori
"""

import time as timing
start = timing.time()

import numpy as np
import pandas as pd

from glob import glob
import netCDF4
import matplotlib.pyplot as plt
import xarray as xr
import gzip

gather = False
#gather = True

campaign = 'tripex'
hydroset = 'all_hydro'
#hydroset = 'only_snow'
#hydroset = 'only_ice'

#hydroset = 'no_snow'
#hydroset = 'only_liquid'
#hydroset = 'only_graupel_hail'

pamtra_radar_data_path = '/data/optimice/pamtra_runs/'+campaign+'/data/'
icon_data_path = '/data/inscape/icon/experiments/juelich/testbed/testbed_'
pluviopath = '/data/data_hatpro/jue/data/pluvio/netcdf/'
mrr_path = '/data/data_hatpro/jue/data/mrr/mrr_ave-6-0-0-6_nc/'

if campaign == 'tripex':
  J10files = sorted(glob(pamtra_radar_data_path + hydroset + '/????????' + hydroset + '_mom_KiXPol.nc'))
  J35files = sorted(glob(pamtra_radar_data_path + hydroset + '/????????' + hydroset + '_mom_Joyrad35.nc'))
  G94files = sorted(glob(pamtra_radar_data_path + hydroset + '/????????' + hydroset + '_mom_Joyrad94.nc'))
elif campaign == 'tripex-pol':
  J10files = sorted(glob(pamtra_radar_data_path + hydroset + '/????????' + hydroset + '_mom_Joyrad10.nc'))
  J35files = sorted(glob(pamtra_radar_data_path + hydroset + '/????????' + hydroset + '_mom_Joyrad35.nc'))
  G94files = sorted(glob(pamtra_radar_data_path + hydroset + '/????????' + hydroset + '_mom_Grarad94.nc'))

available_dates = sorted([i.split('/')[-1][:8] for i in J35files])


if gather:
  for i, date in enumerate(available_dates[:]):
    print(i, date, J10files[i].split('/')[-1])
    iconfile = icon_data_path + date + '/METEOGRAM_patch001_' + date + '_joyce.nc'
    icon_data = netCDF4.Dataset(iconfile).variables
    icon_time = netCDF4.num2date(icon_data['time'][:], icon_data['time'].units)
    icon_numtime = netCDF4.date2num(icon_time,
                                    'seconds since 1970-01-01 00:00:00 UTC')
    runtime = icon_data['time'][:]
    icon_rain = icon_data['RAIN_GSP'][:] + icon_data['RAIN_CON'][:]
    icon_snow = icon_data['SNOW_GSP'][:] + icon_data['RAIN_CON'][:]
    icon_prec = icon_rain + icon_snow
    DF = pd.DataFrame(data=icon_prec, index=icon_time, columns=['ICON_PREC'])
    DF.to_hdf('data/' + campaign + '_' + hydroset + 'precipitation_icon.h5',
              key='stat', mode='a', append=True)
    
    ym = date[2:6]
    pluviofile = pluviopath + ym + '/pluvio2_jue_' + date + '.nc'
    PLUVIO = xr.open_dataset(pluviofile)
    pluvio_data = netCDF4.Dataset(pluviofile).variables
    #pluvio_rrRT = PLUVIO['rain_rate']
    pluvio_accNRT = PLUVIO['r_accum_NRT']*60.0
    pluvio_time = PLUVIO['time']
    nat = PLUVIO['time'].isnull()
    PF = pd.DataFrame(data=pluvio_accNRT[~nat].to_masked_array(),
                      index=pluvio_time[~nat].to_masked_array(),
                      columns=['PLUVIO_PREC'])
    PF.to_hdf('data/' + campaign + '_' + hydroset + 'precipitation_pluvio.h5',
              key='stat', mode='a', append=True)
    
    mrrfile = mrr_path + ym + '/' + date[2:] + '_jue_mrr_ave-6-0-0-6.nc.gz'
    try:
      with gzip.open(mrrfile) as gz:
        mrr = netCDF4.Dataset('dummy', mode='r', memory=gz.read())
    except(FileNotFoundError):
      mrr = netCDF4.Dataset(mrrfile[:-3], mode='r')
    mrr_time = netCDF4.num2date(mrr['time'][:], units='seconds since 1970-1-1')
    mrr_rr = mrr['MRR_RR'][:,0]
    MF = pd.DataFrame(data=mrr_rr, index=mrr_time, columns=['MRR_PREC'])
    MF.to_hdf('data/' + campaign + '_' + hydroset + 'precipitation_mrr.h5',
              key='stat', mode='a', append=True)
    
else:
  print('opening files')
  freq='10min'
  iconfile = 'data/' + campaign + '_' + hydroset + 'precipitation_icon.h5'
  icon = pd.read_hdf(iconfile, key='stat')
  icon = icon.reset_index().drop_duplicates(subset='index',
                                            keep='last').set_index('index')
  pluviofile = 'data/' + campaign + '_' + hydroset + 'precipitation_pluvio.h5'
  pluvio = pd.read_hdf(pluviofile, key='stat')
  mrrfile = 'data/' + campaign + '_' + hydroset + 'precipitation_mrr.h5'
  mrr = pd.read_hdf(mrrfile, key='stat')
  mrr = mrr.resample(freq).apply(np.nanmean)
  
  ixi = pd.date_range(start='2015-11-11', end='2016-1-4', freq='9s')
  icon.reindex(ixi)
  icon = icon.resample(freq).apply(np.nanmean)
  
  ixp = pd.date_range(start='2015-11-11', end='2016-1-4', freq='1min')
  pluvio.reindex(ixp)
  pluvio = pluvio.resample(freq).apply(np.nanmean)
  
  plt.hist(icon['ICON_PREC'], bins=50)
#  plt.figure()
#  plt.hist(pluvio['PLUVIO_PREC'])
#  plt.figure()
  plt.hist(mrr['MRR_PREC'], bins=50)
  

end = timing.time()
print(end-start)
