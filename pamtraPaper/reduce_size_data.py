#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 14:59:21 2019

@author: dori
"""

import xarray as xr
datapath = '/data/inscape/icon/experiments/juelich/testbed/'
filename = 'testbed_20151119/METEOGRAM_patch001_20151119_joyce.nc'
DS = xr.open_dataset(datapath + filename)
delete_vars = ['PEXNER', 'RHO', 'THETAV', 'TKE', 'NIACT', 'QV_DIA', 'QC_DIA',
               'QI_DIA', 'CLC', 'TKVM', 'TKVH', 'PHALF', 'T_SO', 'W_SO',
               'W_SO_ICE', 'PL_COV', 'LA_IND', 'RO_DEPT', 'Z0', 'QV_S', 'W_I',
               'W_SNOW', 'RUNOFF_S', 'RUNOFF_G', 'T_SNOW', 'T_G', 'FRESHSNW',
               'RHO_SNOW', 'H_SNOW', 'FR_SEAICE', 'P_SFC', 'TCM', 'TCH',
               'SHFL', 'LHFL', 'VIO3', 'HMO3', 'T2M', 'TD2M', 'VBMAX10M',
               'dyn_gust', 'con_gust', 'cape_ml', 'SOBT', 'THBT', 'SOBS',
               'THBS', 'ALB', 'RAIN_GSP', 'SNOW_GSP', 'RAIN_CON', 'SNOW_CON',
               'H_ICE', 'CLCT', 'CLCL', 'CLCM', 'CLCH', 'hbas_con', 'htop_con',
               'UMFL_S', 'VMFL_S', 'SWDIFU_S', 'SWDIFD_S', 'PAB_S', 'SWDIR_S',
               'TQV', 'TQC', 'TQI', 'TQR', 'TQS', 'TQV_DIA', 'TQC_DIA',
               'TQI_DIA', 'QV', 'time_step', 'date']
DS.to_netcdf('METEOGRAM_20151119_joyce_full.nc', mode='w', format='NETCDF4',
             engine='netcdf4')
keylist = [i for i in DS.keys()]
DS = DS.drop(delete_vars)
DS.to_netcdf('METEOGRAM_20151119_joyce_reduced.nc', mode='w', format='NETCDF4')
keylist = [i for i in DS.keys()]

DA=DS.isel(time=range(0, len(DS.time), 33)) # 7=63 sec / 33=297sec5min
DA.to_netcdf('METEOGRAM_20151119_joyce_5min.nc', mode='w', format='NETCDF4')

keep = ['time', 'height', 'height_2', 'depth', 'depth_2', 'P', 'T', 'U', 'V',
        'W', 'QC', 'QI', 'QR', 'QS', 'REL_HUM', 'QG', 'QH', 'QNI', 'QNS',
        'QNR', 'QNG', 'QNH', 'QNC', 'T_S', 'U10M', 'V10M']
compression = [{'zlib':True} for i in keep]
DA.to_netcdf('METEOGRAM_20151119_joyce_5min_compressed.nc',
             mode='w', format='NETCDF4',
             encoding=dict(zip(keep, compression)))

DS = DS.drop(keep)
DS.to_netcdf('METEOGRAM_20151119_joyce_none.nc', mode='w', format='NETCDF4')
keylist = [i for i in DS.keys()]
