#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 15:40:31 2019

@author: dori
"""

import time as timer
start = timer.time()

import numpy as np
import pandas as pd

from glob import glob
import netCDF4
import matplotlib.pyplot as plt
#import xarray

campaign = 'tripex'
hydroset = 'all_hydro'

pamtra_radar_data_path = '/data/optimice/pamtra_runs/'+campaign+'/data/'
icon_data_path = '/data/inscape/icon/experiments/juelich/testbed/testbed_'

if campaign == 'tripex':
  J10files = sorted(glob(pamtra_radar_data_path + hydroset + '/????????' + hydroset + '_mom_KiXPol.nc'))
  J35files = sorted(glob(pamtra_radar_data_path + hydroset + '/????????' + hydroset + '_mom_Joyrad35.nc'))
  G94files = sorted(glob(pamtra_radar_data_path + hydroset + '/????????' + hydroset + '_mom_Joyrad94.nc'))
elif campaign == 'tripex-pol':
  J10files = sorted(glob(pamtra_radar_data_path + hydroset + '/????????' + hydroset + '_mom_Joyrad10.nc'))
  J35files = sorted(glob(pamtra_radar_data_path + hydroset + '/????????' + hydroset + '_mom_Joyrad35.nc'))
  G94files = sorted(glob(pamtra_radar_data_path + hydroset + '/????????' + hydroset + '_mom_Grarad94.nc'))

available_dates = sorted([i.split('/')[-1][:8] for i in J10files])

def running_mean(x, N, minN):
    csumnan = np.cumsum((~np.isfinite(np.insert(x, 0, 0))).astype(int))
    nannum = csumnan[N:] - csumnan[:-N]
    mask = ((N-nannum) > minN)
    Filter = mask.astype(float)
    Filter[~mask] = np.nan
    cumsum = np.nancumsum(np.insert(x, 0, 0)) 
    return Filter * (cumsum[N:] - cumsum[:-N]) / (float(N)-nannum)

def running_mean_2d(xx, N, minN=0):
  add = int((N - (N & 1))/2)
  addition = np.zeros([xx.shape[0], add])*np.nan
  x = np.concatenate([addition, xx, addition], axis=1)
  csumnan = np.cumsum((~np.isfinite(np.insert(x, 0, 0, axis=1))).astype(int),
                      axis=1)
  nannum = csumnan[:, N:] - csumnan[:, :-N]
  mask = ((N-nannum) >= minN)
  Filter = mask.astype(float)
  Filter[~mask] = np.nan
  csum = np.nancumsum(np.insert(x, 0, 0, axis=1), axis=1)
  return Filter * (csum[:, N:] - csum[:, :-N]) / (float(N)-nannum)

if __name__ == '__main__':
  for i, date in enumerate(available_dates[:]):
    print(i, date)
    #radar_file = '/data/optimice/tripex/tripex_level_02_samd/tripex_joy_tricr00_l2_any_v00_' + date + '000000.nc'
    radar_file = '/data/optimice/tripex/tripex_level_02_test/tripex_joy_tricr00_l2_any_v00_' + date + '000000.nc'
    Nt = netCDF4.Dataset(radar_file).dimensions.get('time').size
    Nr = netCDF4.Dataset(radar_file).dimensions.get('range').size
    print(radar_file)
    radar_data = netCDF4.Dataset(radar_file).variables
    radar_time = netCDF4.num2date(radar_data['time'][:], radar_data['time'].units)

    DF = pd.DataFrame()
       
    DF['Z10'] = radar_data['dbz_x'][:].T.flatten()
    DF['V10'] = radar_data['rv_x'][:].T.flatten()
    DF['W10'] = radar_data['sw_x'][:].T.flatten()
    DF['Z35'] = radar_data['dbz_ka'][:].T.flatten()
    DF['V35'] = radar_data['rv_ka'][:].T.flatten()
    DF['W35'] = radar_data['sw_ka'][:].T.flatten()
    DF['Z94'] = radar_data['dbz_w'][:].T.flatten()
    DF['V94'] = radar_data['rv_w'][:].T.flatten()
    DF['W94'] = radar_data['sw_w'][:].T.flatten()
    
    V10=radar_data['rv_x'][:].T
    V35=radar_data['rv_ka'][:].T
    V94=radar_data['rv_w'][:].T
    DF['V10avg']=running_mean_2d(V10, 299, 75).flatten() # 299 * 4 sec = 20 min
    DF['V35avg']=running_mean_2d(V35, 299, 75).flatten() # 75 * 4 sec = 5 min of measurements
    DF['V94avg']=running_mean_2d(V94, 299, 75).flatten()
    
    time = radar_data['time'][:] - netCDF4.date2num(pd.to_datetime(date),
                                                    'seconds since 1970-01-01 00:00:00 UTC')
    DF['runtime'] = np.tile(time,[Nr,1]).flatten()
    DF['Hgt'] = np.tile(radar_data['range'][:][np.newaxis].T,[1,Nt]).flatten() + 112.5
    DF['P'] = radar_data['Pres_Cl'][:].T.flatten()
    DF['T'] = radar_data['Temp_Cl'][:].T.flatten()
    DF['RH'] = radar_data['RelHum_Cl'][:].T.flatten()

    DF['quality_x'] = radar_data['quality_flag_offset_x'][:].T.flatten()
    DF['quality_w'] = radar_data['quality_flag_offset_w'][:].T.flatten()
    
    DF.dropna(how='all', subset=['Z10', 'Z35', 'Z94'], inplace=True)
    DF.to_hdf('data/radar/' + campaign + '_data_radar_avg.h5',
              key='stat', mode='a', append=True)
    for col in DF.columns:
      DF[col].to_hdf('data/radar/' + campaign + '_' + col + '_data_radar.h5',
                     key='stat', mode='a', append=True)
    #xr = xarray.Dataset.from_dataframe(DF)
    #key_list = [i for i in xr.keys()]
    #compress = [{'zlib':True} for i in key_list]
    #xr.to_netcdf('data/radar/radar_compress.nc', mode='w', format='NETCDF4',
    #             encoding=dict(zip(key_list, compress)))
    #xr.to_netcdf('data/radar/radar.nc', mode='w', format='NETCDF4')
  #end = timer.time()
  #print(end-start)
  #DF2 = pd.read_hdf('data/' + campaign + '_data_radar_avg.h5', key='stat')
  end = timer.time()
  print(end-start)