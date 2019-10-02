#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 16:54:03 2019

@author: dori
"""

import pandas as pd
import numpy as np

T0 = 273.15

radarvars = ['Z10', 'Z35', 'Z94', 'N10', 'N35', 'N94', 'S10', 'S35', 'S94',
             'V10', 'V35', 'V94', 'W10', 'W35', 'W94',
             'V10avg', 'V35svg', 'V94svg']

def read_variables(path='data/', hydroset='all_hydro', suffix='pamtra_icon.h5',
                   varlist=['runtime'], pamtra=False,
                   attATM=False, attHYD=False,
                   minhour=0.0, maxhour=9e9):
  series_list = []
  if 'runtime' not in varlist:
    varlist.append('runtime')
  for var in varlist:
    hydro = hydroset
    if len(hydroset):
      hydro = hydroset + '_'
    filename = path + 'tripex_' + hydro + var + '_data_' + suffix
    series_list.append(pd.read_hdf(filename, key='stat'))
  radarvarlist = list(set(radarvars) & set(varlist))
  data = pd.concat(series_list, axis=1)
  if radarvarlist:
    data[radarvarlist] = data[radarvarlist].replace(-9999.0, np.nan)
  if pamtra:
    if 'T' in varlist:
      data['T'] = data['T'] - 273.15
    if 'V10' in varlist:
      data['V10'] = -1.0*data['V10']
    if 'V35' in varlist:
      data['V35'] = -1.0*data['V35']
    if 'V94' in varlist:
      data['V94'] = -1.0*data['V94']
  
  if pamtra:
    if attATM:
      if ('Z10' in varlist):
        data['Z10'] = data['Z10'] - data['A10']
      if ('Z35' in varlist):
        data['Z35'] = data['Z35'] - data['A35']
      if ('Z94' in varlist):
        data['Z94'] = data['Z94'] - data['A94']
    if attHYD:
      if ('Z10' in varlist):
        data['Z10'] = data['Z10'] - data['H10']
      if ('Z35' in varlist):
        data['Z35'] = data['Z35'] - data['H35']
      if ('Z94' in varlist):
        data['Z94'] = data['Z94'] - data['H94']
      

  if ('Z10' in varlist) and ('Z35' in varlist):
    data['DWRxk'] = data['Z10'] - data['Z35']
  if ('Z94' in varlist) and ('Z35' in varlist):
    data['DWRkw'] = data['Z35'] - data['Z94']
  return slice_data(data, 'runtime', 3600.*minhour, 3600.*maxhour)


def read_prepare(hydroset='all_hydro', minhour=0.0, maxhour=999999.9, suffix='', filename=None):
  if filename is None:
    filename = 'data/'+hydroset+'data_pamtra_icon'+suffix+'.h5'
  data = pd.read_hdf(filename, key='stat')
  data['T'] = data['T'] - 273.15
  data['V10'] = -1.0*data['V10']
  data['V35'] = -1.0*data['V35']
  data['V94'] = -1.0*data['V94']
  data[data==-9999] = np.nan
  #data.dropna(inplace=True)
  data['DWRxk'] = data['Z10'] - data['Z35']
  data['DWRkw'] = data['Z35'] - data['Z94']
  #data = data[(data['runtime']>3600.*minhour)*(data['runtime']<3600.*maxhour)]
  return slice_data(data, 'runtime', 3600.*minhour, 3600.*maxhour)


def read_radar(campaign, cols=[], minhour=0.0, maxhour=999999.9, avg='', filename=None):
  if 'runtime' not in cols:
    cols.append('runtime')
  if len(cols)==1:
    cols=None
  if filename is None:
    filename = 'data/'+campaign+'_data_radar'+avg+'.h5'
  data = pd.read_hdf(filename, key='stat', columns=cols)
  #data.dropna(inplace=True)
  cols=data.columns
  if ('Z10' in cols) and ('Z35' in cols):
    data['DWRxk'] = data['Z10'] - data['Z35']
  if ('Z94' in cols) and ('Z35' in cols):
    data['DWRkw'] = data['Z35'] - data['Z94']
  #data = data[(data['runtime']>3600.*minhour)*(data['runtime']<3600.*maxhour)]
  return slice_data(data, 'runtime', 3600.*minhour, 3600.*maxhour)


def slice_data(data, varname, minvalue=-np.inf, maxvalue=np.inf,
               left=False, right=False):
  minmask = (data[varname]>=minvalue) if left else (data[varname]>minvalue)
  maxmask = (data[varname]<=maxvalue) if right else (data[varname]<maxvalue)
  return data[minmask&maxmask]


icon150heights = [21000, 20617.7923391807, 20290.2350703975, 19982.3469255595, 
    19687.435072837, 19402.3933650287, 19125.4122583281, 18855.3018784347, 
    18591.2174854262, 18332.526574574, 18078.7366066168, 17829.4523000272, 
    17584.3487582748, 17343.153708233, 17105.6352855012, 16871.5933550092, 
    16640.8531719795, 16413.2606425969, 16188.6787084653, 15967.5088055511, 
    15748.6404443982, 15532.4519150312, 15318.8516703129, 15107.7558114719, 
    14899.08712008, 14692.7742508499, 14488.7510529727, 14286.9559952076, 
    14087.3316754788, 13889.8243998731, 13694.3838190755, 13500.9626126768, 
    13309.5162136479, 13120.0025667213, 12932.3819155615, 12746.6166145103, 
    12562.6709614159, 12380.5110486355, 12200.1046297755, 12021.4210001143, 
    11844.4308889707, 11669.1063625404, 11495.420735938, 11323.348493363, 
    11152.8652154582, 10983.947513054, 10816.5729666031, 10650.7200706956, 
    10486.3681831264, 10323.4974780493, 10162.0889028117, 10002.1241381104, 
    9843.58556115232, 9686.45621154006, 9530.71975963426, 9376.36047717196, 
    9223.36320994482, 9071.71335236187, 8921.39682374051, 8772.40004618573, 
    8624.70992393207, 8478.31382403597, 8333.19955831737, 8189.35536645983, 
    8046.76990018764, 7905.4322084464, 7765.33172352109, 7626.45824803236, 
    7488.80194275746, 7352.35331522815, 7217.10320906248, 7083.04279399218, 
    6950.16355655133, 6818.45729139588, 6687.91609322729, 6558.53234929653, 
    6430.29873246809, 6303.20819482616, 6177.25396180829, 6052.42952685405, 
    5928.72864655888, 5806.14533632576, 5684.67386650935, 5564.30875905009, 
    5445.04478459738, 5326.87696012367, 5209.80054703349, 5093.81104977361, 
    4978.90421495311, 4865.07603098453, 4752.32272826002, 4640.64077987891, 
    4530.0269029464, 4420.478060466, 4311.99146385187, 4204.56457609104, 
    4098.19511558936, 3992.88106073987, 3888.62065525708, 3785.41241432627, 
    3683.25513162343, 3582.14788726808, 3482.09005677969, 3383.08132111695, 
    3285.12167788982, 3188.21145384588, 3092.35131874567, 2997.54230075775, 
    2903.78580352154, 2811.08362504686, 2719.43797864326, 2628.85151610001, 
    2539.32735337076, 2450.86909905519, 2363.48088601599, 2277.16740652332, 
    2191.9339513835, 2107.78645358588, 2024.7315370947, 1942.77657152497, 
    1861.9297335781, 1782.20007627993, 1703.59760726952, 1626.13337764181, 
    1549.81958316532, 1474.66968009583, 1400.69851831368, 1327.92249516096, 
    1256.35973419083, 1186.03029413285, 1116.95641481621, 1049.16280871212, 
    982.677009351001, 917.529791428042, 853.755682366588, 791.393592134995, 
    730.487598267305, 671.087938038428, 613.252282452167, 557.047402048352, 
    502.551391351508, 449.856713567836, 399.074492338514, 350.340780758514, 
    303.82613180232, 259.751053575898, 218.412893285042, 180.237701543821, 
    145.897229870517, 116.815557244347, 96.8155572443466]