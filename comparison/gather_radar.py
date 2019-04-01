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
    DF.to_hdf('data/' + campaign + '_data_radar.h5',
              key='stat',mode='a',append=True)
    
  end = timer.time()
  print(end-start)