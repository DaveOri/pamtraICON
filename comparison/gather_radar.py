#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 15:40:31 2019

@author: dori
"""

import time as timer
import numpy as np
import pandas as pd
from glob import glob
import netCDF4
import matplotlib.pyplot as plt
from READ import icon150heights

start = timer.time()
icon150heights = icon150heights[::-1]
Hmax = np.max(np.array(icon150heights))
Hmin = np.min(np.array(icon150heights))
heightCenter = 0.5*(np.array(icon150heights)[1:]+np.array(icon150heights)[:-1])
Tmin = 1
Tmax = 86400
deltat = 9
iconTimes = np.arange(4.5, 86405, deltat)
timeCenters = np.arange(0+deltat, Tmax+deltat, deltat)


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


def find_bin(x, bins, binCenter):
  idx = np.digitize(x, bins)
  return binCenter[idx-1]


def find_height_Icon(x):
  return find_bin(x, np.array(icon150heights), heightCenter)


def find_time_Icon(x):
  return find_bin(x, iconTimes, timeCenters)


def dB(x):
  return 10.0*np.log10(x)


def Bd(x):
  return 10.0**(0.1*x)


def linearMean(x):
  return dB(np.nanmean(Bd(x)))


def linearWeightedMean(x, weights):
  w = Bd(x)
  return np.nansum(x*w)/np.nansum(w)


aggDict = {'Z10':linearMean,
           'V10':np.nanmean,
           'W10':np.nanmean,
           'Z35':linearMean,
           'V35':np.nanmean,
           'W35':np.nanmean,
           'Z94':linearMean,
           'V94':np.nanmean,
           'W94':np.nanmean,
           'V10avg':np.nanmean,
           'V35avg':np.nanmean,
           'V94avg':np.nanmean,
           'runtime':np.nanmean,
           'unixtime':np.nanmean,
           'Hgt':np.nanmean,
           'P':np.nanmean,
           'T':np.nanmean,
           'RH':np.nanmean,
           'quality_x':np.max,
           'quality_w':np.max,
           'groupH':np.nanmean}


def reduction(x):
  d = {}
  d['Z10'] = linearMean(x['Z10'])
  d['V10'] = linearWeightedMean(x['V10'], x['Z10'])
  d['W10'] = linearWeightedMean(x['W10'], x['Z10'])
  d['Z35'] = linearMean(x['Z35'])
  d['V35'] = linearWeightedMean(x['V35'], x['Z35'])
  d['W35'] = linearWeightedMean(x['W35'], x['Z35'])
  d['Z94'] = linearMean(x['Z94'])
  d['V94'] = linearWeightedMean(x['V94'], x['Z94'])
  d['W94'] = linearWeightedMean(x['W94'], x['Z94'])
  d['V10avg'] = linearWeightedMean(x['V10avg'], x['Z10'])
  d['V35avg'] = linearWeightedMean(x['V35avg'], x['Z35'])
  d['V94avg'] = linearWeightedMean(x['V94avg'], x['Z94'])
  d['runtime'] = np.nanmean(x['runtime'])
  d['unixtime'] = np.nanmean(x['unixtime'])
  d['Hgt'] = np.nanmean(x['Hgt'])
  d['P'] = np.nanmean(x['P'])
  d['T'] = np.nanmean(x['T'])
  d['RH'] = np.nanmean(x['RH'])
  d['quality_x'] = np.max(x['quality_x'])
  d['quality_w'] = np.max(x['quality_w'])
  return pd.Series(d, index=d.keys())


def running_mean(x, N, minN):
    csumnan = np.cumsum((~np.isfinite(np.insert(x, 0, 0))).astype(int))
    nannum = csumnan[N:] - csumnan[:-N]
    mask = ((N-nannum) > minN)
    Filter = mask.astype(float)
    Filter[~mask] = np.nan
    cumsum = np.nancumsum(np.insert(x, 0, 0)) 
    return Filter * (cumsum[N:] - cumsum[:-N]) / (float(N)-nannum)


def running_mean_2d(xx, N, minN=0, weights=None):
  if weights is None:
    weights = np.ones(xx.shape)
  xx = xx*weights
  add = int((N - (N & 1))/2)
  addition = np.zeros([xx.shape[0], add])*np.nan
  x = np.concatenate([addition, xx, addition], axis=1)
  w = np.concatenate([addition, weights, addition], axis=1)
  csumnan = np.cumsum((~np.isfinite(np.insert(x, 0, 0, axis=1))).astype(int),
                      axis=1)
  nannum = csumnan[:, N:] - csumnan[:, :-N]
  mask = ((N-nannum) >= minN)
  Filter = mask.astype(float)
  Filter[~mask] = np.nan
  csum = np.nancumsum(np.insert(x, 0, 0, axis=1), axis=1)
  wsum = np.nancumsum(np.insert(w, 0, 0, axis=1), axis=1)
  return Filter * (csum[:, N:] - csum[:, :-N]) / (wsum[:, N:] - wsum[:, :-N])
  #return Filter * (csum[:, N:] - csum[:, :-N]) / (float(N)-nannum)
  

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
    Z10 = radar_data['dbz_x'][:].T
    Z35 = radar_data['dbz_ka'][:].T
    Z94 = radar_data['dbz_w'][:].T
    DF['Z10'] = Z10.flatten()
    DF['V10'] = radar_data['rv_x'][:].T.flatten()
    DF['W10'] = radar_data['sw_x'][:].T.flatten()
    DF['Z35'] = Z35.flatten()
    DF['V35'] = radar_data['rv_ka'][:].T.flatten()
    DF['W35'] = radar_data['sw_ka'][:].T.flatten()
    DF['Z94'] = Z94.flatten()
    DF['V94'] = radar_data['rv_w'][:].T.flatten()
    DF['W94'] = radar_data['sw_w'][:].T.flatten()
    
    V10=radar_data['rv_x'][:].T
    V35=radar_data['rv_ka'][:].T
    V94=radar_data['rv_w'][:].T
    DF['V10avg']=running_mean_2d(V10, 299, 75, Bd(Z10)).flatten() # 299 * 4 sec = 20 min
    DF['V35avg']=running_mean_2d(V35, 299, 75, Bd(Z35)).flatten() # 75 * 4 sec = 5 min of measurements
    DF['V94avg']=running_mean_2d(V94, 299, 75, Bd(Z94)).flatten()
    
    time = radar_data['time'][:] - netCDF4.date2num(pd.to_datetime(date),
                                                    'seconds since 1970-01-01 00:00:00 UTC')
    DF['runtime'] = np.tile(time,[Nr,1]).flatten()
    DF['unixtime'] = np.tile(radar_data['time'][:],[Nr,1]).flatten()
    DF['Hgt'] = np.tile(radar_data['range'][:][np.newaxis].T,[1,Nt]).flatten() + 112.5
    DF['P'] = radar_data['Pres_Cl'][:].T.flatten()
    DF['T'] = radar_data['Temp_Cl'][:].T.flatten()
    DF['RH'] = radar_data['RelHum_Cl'][:].T.flatten()

    DF['quality_x'] = radar_data['quality_flag_offset_x'][:].T.flatten()
    DF['quality_w'] = radar_data['quality_flag_offset_w'][:].T.flatten()
    
    mask = (DF['Hgt']>=Hmax)+(DF['Hgt']<=Hmin)
    mask = mask + (DF['runtime']>=Tmax)+(DF['runtime']<=Tmin)
    DF.loc[mask, ['Z10','Z35','Z94']] = np.nan
    DF.dropna(how='all', subset=['Z10', 'Z35', 'Z94'], inplace=True)
    
    #
    DF['groupT'] = find_time_Icon(DF.loc[:,'runtime'])
    DF['groupH'] = find_height_Icon(DF.loc[:,'Hgt'])
    groups = DF.groupby(['groupH', 'groupT'])
    # rDF = groups.aggregate(aggDict)
    rDF = groups.apply(reduction)
    rDF.dropna(how='all', subset=['Z10','Z35','Z94'])
    # rDF.drop('groupH', axis=1, inplace=True)
    #
    
    DF.to_hdf('data/radar/' + campaign + '_data_radar_avg.h5',
              key='stat', mode='a', append=True)
    rDF.to_hdf('data/radar/' + campaign + '_data_radar_avg_regrid.h5',
               key='stat', mode='a', append=True)
    for col in rDF.columns:
      DF[col].to_hdf('data/radar/' + campaign + '_' + col + '_data_radar.h5',
                     key='stat', mode='a', append=True)
      rDF[col].to_hdf('data/radar/' + campaign + '_' + col + '_data_radar_regrid.h5',
                      key='stat', mode='a', append=True)

      
  end = timer.time()
  print(end-start)
  