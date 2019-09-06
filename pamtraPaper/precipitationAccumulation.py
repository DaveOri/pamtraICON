#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 13:23:57 2019

@author: dori
"""

import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')

datapath = '/data/inscape/icon/experiments/juelich/testbed/'
filename = 'testbed_20151119/METEOGRAM_patch001_20151119_joyce.nc'
ICON = xr.open_dataset(datapath + filename)
icon_rain = ICON['RAIN_GSP'] + ICON['RAIN_CON']
icon_snow = ICON['SNOW_GSP'] + ICON['RAIN_CON']

pluviopath = '/data/data_hatpro/jue/data/pluvio/netcdf/'
pluviofile = pluviopath + '1511/pluvio2_jue_20151119.nc'
PLUVIO = xr.open_dataset(pluviofile)
pluvio_TNRT = PLUVIO['r_accum_NRT']

plt.plot(ICON['time'], icon_rain, label='rain')
plt.plot(ICON['time'], icon_snow, label='snow')
plt.plot(ICON['time'], icon_snow+icon_rain, label='total icon')
nat = PLUVIO['time'].isnull()
plt.plot(PLUVIO['time'][~nat],
         (pluvio_TNRT[~nat]).cumsum(),
         label='pluvio')
plt.legend()