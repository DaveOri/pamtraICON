#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 11:48:30 2019

@author: dori
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import sys
sys.path.append('..')
from READ import slice_data
from READ import read_variables
from statistic import hist_and_plot

varlist = ['Hgt', 'P', 'T', 'RH', 'Z10', 'Z35', 'Z94', 'N10', 'N35', 'N94',
           'QI', 'QC', 'QG', 'QH', 'QR', 'QS',
           'QNI', 'QNC', 'QNG', 'QNH', 'QNR', 'QNS']

hydro = ['QS', 'QI', 'QC', 'QG', 'QH', 'QR',
         'QNS', 'QNI', 'QNC', 'QNG', 'QNH', 'QNR']


pamtra = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='all_hydro', suffix='pamtra_icon.h5',
                        varlist=varlist, minhour=6.0)

data = pamtra.dropna(subset=['DWRxk'])
data['RI'] = data['QI']/data['QNI']
data['RS'] = data['QS']/data['QNS']
data['RC'] = data['QC']/data['QNC']
data['RR'] = data['QR']/data['QNR']
data['RG'] = data['QG']/data['QNG']
data['RH'] = data['QH']/data['QNH']

data[hydro] = np.log10(data[hydro])

data['QNS'][data['QNS']==-np.inf] = -302
data['QNI'][data['QNI']==-np.inf] = -289
data['QNR'][data['QNR']==-np.inf] = -219
data['QNC'][data['QNC']==-np.inf] = -2.5
data['QNH'][data['QNH']==-np.inf] = -304
data['QNG'][data['QNG']==-np.inf] = -304
data['QS'][data['QS']==-np.inf] = -307
data['QI'][data['QI']==-np.inf] = -301
data['QR'][data['QR']==-np.inf] = -224
data['QC'][data['QC']==-np.inf] = -10
data['QH'][data['QH']==-np.inf] = -307
data['QG'][data['QG']==-np.inf] = -307


def histogram_data_dwr_2(data, col):
  ax = data.hist(column=col, bins=100, density=True, alpha=0.5)
  data[data['DWRxk']<-2].hist(column=col, ax=ax, bins=100, density=True, alpha=0.5)
  ax = plt.gca()
  ax.legend(labels=['all data', 'DWRxk<-2'])

for var in varlist:
  histogram_data_dwr_2(data, var)

data['RS'][data['RS']==np.inf] = 0.000021
data['RI'][data['RI']==np.inf] = 0.0000011
data['RR'][data['RR']==np.inf] = 0.0000032
data['RC'][data['RC']==np.inf] = np.nan#0.0000066
data['RH'][data['RH']==np.inf] = 0.0006
data['RG'][data['RG']==np.inf] = 0.0006
histogram_data_dwr_2(data, 'RI')
histogram_data_dwr_2(data, 'RS')
histogram_data_dwr_2(data, 'RC')
histogram_data_dwr_2(data, 'RR')
histogram_data_dwr_2(data, 'RG')
histogram_data_dwr_2(data, 'RH')
