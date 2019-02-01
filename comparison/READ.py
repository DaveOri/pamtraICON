#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 16:54:03 2019

@author: dori
"""

import pandas as pd
import numpy as np

T0 = 273.15

def read_prepare(hydroset='all_hydro', minhour=0.0, maxhour=999999.9 ):
  data = pd.read_hdf('data/'+hydroset+'data_pamtra_icon.h5',key='stat')
  data['T'] = data['T'] - 273.15
  data['V10'] = -1.0*data['V10']
  data['V35'] = -1.0*data['V35']
  data['V94'] = -1.0*data['V94']
  data[data==-9999] = np.nan
  data.dropna(inplace=True)
  data['DWRxk'] = data['Z10'] - data['Z35']
  data['DWRkw'] = data['Z35'] - data['Z94']
  data = data[(data['runtime']>3600.*minhour)*(data['runtime']<3600.*maxhour)]
  return data