#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 16:00:12 2019

@author: dori
"""

import matplotlib.pyplot as plt
plt.close('all')
import matplotlib.dates as md
import pandas as pd
import numpy as np
import netCDF4 as nc

pamX = nc.Dataset('all_hydro/20151119all_hydro_spe_KiXPol.nc')
pamK = nc.Dataset('all_hydro/20151119all_hydro_spe_Joyrad35.nc')
pamW = nc.Dataset('all_hydro/20151119all_hydro_spe_Joyrad94.nc')
pamP = nc.Dataset('20151119hatpro.nc')

rad3 = nc.Dataset('/data/optimice/tripex/tripex_level_02_test/tripex_joy_tricr00_l2_any_v00_20151119000000.nc')
hatp = nc.Dataset('/data/hatpro/jue/data/level1/1511/sups_joy_mwr00_l1_tb_v01_20151119000003.nc')

xfmt = md.DateFormatter('%m-%d %H')
xfmt = md.DateFormatter('%H')
ylim=(0,12000)
xDataLim = -1
figsize31=(18,18)
figsize21=(18,12)
versus = -1 # Top Down
versus =  1 # Bottom Up

# Extract netCDF4 time variables
def getTime(ncdata, timevar):
  return nc.num2date(times=ncdata.variables[timevar][:],
                     units=ncdata.variables[timevar].units).squeeze()

# Define Plotting Function
def plot_variable(x,y,v,axes,
                  xlab=None,ylab=None,vlab=None,title=None,
                  vmin=None,vmax=None,xlim=None,ylim=None,
                  cmap='jet'):
    mesh = axes.pcolormesh(x,y,v,vmin=vmin,vmax=vmax,cmap=cmap)
    if title is not None:
        axes.text(0.1,0.9,title,transform=axes.transAxes,weight='black',
                  bbox=dict(facecolor='white'))
    plt.colorbar(mesh,label=vlab,ax=axes)
    if xlab is not None:
        axes.set_xlabel(xlab)
    if ylab is not None:
        axes.set_ylabel(ylab)
    axes.set_xlim(xlim)
    axes.set_ylim(ylim)
  
Hx, ttx, Ax, Zex, MDVx, SWx = readPamtra_nc(pamX)
Ha, tta, Aa, Zea, MDVa, SWa = readPamtra_nc(pamK)
Hw, ttw, Aw, Zew, MDVw, SWw = readPamtra_nc(pamW)

hlim = [0, 10]
TlimK = [0, 150]
TlimV = [100, 300]
vlim = [-8, 8]
splim = [-90, 10]
Zmin, Zmax = -35, 25
Vmin, Vmax = -1, 5
Dmin, Dmax = -5, 20

tidx = 939#6572#4800
hidx = 29#70
selTime = tta[tidx, hidx]
selTS = pd.to_datetime(selTime)
selHeight = Ha[tidx, hidx]

rad94file = '/data/hatpro/jue/data/joyrad94/l0/' + str(selTS.year) \
            + str(selTS.month) + '/' + str(selTS.day) + '/joyrad94_joyce_' \
            + str(selTS.year) + str(selTS.month) + str(selTS.day) \
            + str(selTS.hour) + '.nc'
radar94 = nc.Dataset(rad94file)
rad94var = radar94.variables
t94 = nc.num2date(rad94var['time'][:], 'seconds since 2001-01-01')
t94idx = np.argmin(np.abs(t94 - selTS))
h94idx = np.argmin(np.abs(rad94var['range'][:] - selHeight))
spec94 = rad94var['spec'][t94idx,:,:]
spec94[spec94 == -999] = np.nan
chirp_idx = list(rad94var['range_offsets'][:] - 1)
chirp_idx.append(spec94.shape[0])
v94 = np.zeros(spec94.shape)
for i,j in enumerate(chirp_idx[:-1]):
  v94[j:chirp_idx[i+1],:] = np.tile(rad94var['velocity'][i,:], [chirp_idx[i+1]-j,1])
r94 = np.tile(rad94var['range'][:][:, np.newaxis], spec94.shape[1])
v94 = -v94

rad35file = '/net/ora/20151119mira36spectra.nc'
radar35 = nc.Dataset(rad35file)
rad35var = radar35.variables
t35 = nc.num2date(rad35var['time'][:], 'seconds since 1970-01-01 00:00 UTC')
t35idx = 50#np.argmin(np.abs(t35 - selTS))
r35 = rad35var['range'][:]
h35idx = np.argmin(np.abs(r35 - selHeight))

SNRCorFaCo = rad35var['SNRCorFaCo'][t35idx,:]
radConst = rad35var['RadarConst'][t35idx]
npw = rad35var['npw1'][t35idx]
cal = radConst*SNRCorFaCo*(r35/5000.)**2/npw
spec35 = 10.*np.log10(rad35var['SPCco'][t35idx,:,:]*cal[:, np.newaxis])
v35 = sorted(rad35var['doppler'][:])