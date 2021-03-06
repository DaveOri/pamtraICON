#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 15:40:31 2019

@author: dori
"""

import time as timing
start = timing.time()

import numpy as np
import pandas as pd

from glob import glob
import netCDF4
import matplotlib.pyplot as plt

campaign = 'tripex'

hydroset = 'all_hydro'
hydroset = 'only_snow'
hydroset = 'only_ice'
hydroset = 'no_snow'
hydroset = 'only_liquid'
hydroset = 'only_graupel_hail'

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

available_dates = sorted([i.split('/')[-1][:8] for i in J35files])

def swap_gridx_datatime(dataarray):
  dt = dataarray.rename({'grid_x':'datetime'})
  dt.coords['datetime'].values = dt.datatime.values[:,0]
  return dt

def get_moments(R):
  Z = R['Ze'][:,0,::-1,0,0,0].T
  V = R['Radar_MeanDopplerVel'][:,0,::-1,0,0,0].T
  W = R['Radar_SpectrumWidth'][:,0,::-1,0,0,0].T
  S = R['Radar_Skewness'][:,0,::-1,0,0,0].T
  N = R['Radar_SNR'][:,0,::-1,0,0,0].T
  aH = 2.0*R['Attenuation_Hydrometeors'][:,0,:,0,0].T.cumsum(axis=1)
  aA = 2.0*R['Attenuation_Atmosphere'][:,0,:,0].T.cumsum(axis=1)
  return Z, V, W, S, N, aH, aA

if __name__ == '__main__':

  for i, date in enumerate(available_dates[:]):
    print(i, date, J10files[i].split('/')[-1],
          J35files[i].split('/')[-1], G94files[i].split('/')[-1])
    J10 = netCDF4.Dataset(J10files[i], mode='r').variables
    J35 = netCDF4.Dataset(J35files[i], mode='r').variables
    G94 = netCDF4.Dataset(G94files[i], mode='r').variables
    Z10, V10, W10, S10, N10, H10, A10 = get_moments(J10)
    Z35, V35, W35, S35, N35, H35, A35 = get_moments(J35)
    Z94, V94, W94, S94, N94, H94, A94 = get_moments(G94)
    
    iconfile = icon_data_path + date + '/METEOGRAM_patch001_' + date + '_joyce.nc'
    icon_data = netCDF4.Dataset(iconfile).variables
    icon_time = netCDF4.num2date(icon_data['time'][:], icon_data['time'].units)
    icon_time = netCDF4.date2num(icon_time, 'seconds since 1970-01-01 00:00:00 UTC')
    icon_heights = icon_data['height_2'][:]
    P = icon_data['P'][:].T
    T = icon_data['T'][:].T
    RH = icon_data['REL_HUM'][:].T
    Hgt = np.tile(np.tile(icon_heights.T,[1,1]).T,[1,P.shape[1]])
    time = np.tile(icon_time, [P.shape[0],1])
    runtime = np.tile(icon_data['time'][:],[P.shape[0],1])
    
    
    DF = pd.DataFrame()
    DF['runtime'] = runtime.flatten()
    DF['unixtime'] = time.flatten()
    DF['Z10'] = Z10.flatten()
    DF['V10'] = V10.flatten()
    DF['W10'] = W10.flatten()
    DF['S10'] = S10.flatten()
    DF['N10'] = N10.flatten()
    DF['H10'] = H10.flatten()
    DF['A10'] = A10.flatten()
    
    DF['Z35'] = Z35.flatten()
    DF['V35'] = V35.flatten()
    DF['W35'] = W35.flatten()
    DF['S35'] = S35.flatten()
    DF['N35'] = N35.flatten()
    DF['H35'] = H35.flatten()
    DF['A35'] = A35.flatten()
    
    DF['Z94'] = Z94.flatten()
    DF['V94'] = V94.flatten()
    DF['W94'] = W94.flatten()
    DF['S94'] = S94.flatten()
    DF['N94'] = N94.flatten()
    DF['H94'] = H94.flatten()
    DF['A94'] = A94.flatten()
    
    DF['Hgt'] = Hgt.flatten()
    DF['P'] = P.flatten()
    DF['T'] = T.flatten()
    DF['RH'] = RH.flatten()
    
    DF['QNI'] = icon_data['QNI'][:].T.flatten()
    DF['QNS'] = icon_data['QNS'][:].T.flatten()
    DF['QNG'] = icon_data['QNG'][:].T.flatten()
    DF['QNH'] = icon_data['QNH'][:].T.flatten()
    DF['QNC'] = icon_data['QNC'][:].T.flatten()
    DF['QNR'] = icon_data['QNR'][:].T.flatten()
    DF['QI'] = icon_data['QI'][:].T.flatten()
    DF['QS'] = icon_data['QS'][:].T.flatten()
    DF['QG'] = icon_data['QG'][:].T.flatten()
    DF['QH'] = icon_data['QH'][:].T.flatten()
    DF['QC'] = icon_data['QC'][:].T.flatten()
    DF['QR'] = icon_data['QR'][:].T.flatten()
    
    DF.to_hdf('data/pamtra/' + campaign + '_' + hydroset + '_data_pamtra_icon.h5',
              key='stat',mode='a',append=True)
    for col in DF.columns:
      DF[col].to_hdf('data/pamtra/' + campaign + '_' + hydroset + '_' + col + '_data_pamtra_icon.h5',
                     key='stat',mode='a',append=True)
    
  # f, (ax1,ax2,ax3) = plt.subplots(3, 1, sharex=True, figsize=(18, 18))
  # mesh1 = ax1.pcolormesh(icon_time, icon_heights, P)
  # mesh2 = ax2.pcolormesh(icon_time, icon_heights, T)
  # mesh3 = ax3.pcolormesh(icon_time, icon_heights, RH)
  # plt.colorbar(mesh1, label='P   [Pa]', ax=ax1)
  # plt.colorbar(mesh2, label='T   [K]', ax=ax2)
  # plt.colorbar(mesh3, label='RH   [%]', ax=ax3)

  # f, (ax1,ax2,ax3) = plt.subplots(3, 1, sharex=True, figsize=(18, 18))
  # mesh1 = ax1.pcolormesh(icon_time, icon_heights, Z10)
  # mesh2 = ax2.pcolormesh(icon_time, icon_heights, V10)
  # mesh3 = ax3.pcolormesh(icon_time, icon_heights, S10)
  # plt.colorbar(mesh1, label='Z    [dBZ]', ax=ax1)
  # plt.colorbar(mesh2, label='V    [m/s]', ax=ax2)
  # plt.colorbar(mesh3, label='SW   [m/s]', ax=ax3)

  end = timing.time()
  print(end-start)
