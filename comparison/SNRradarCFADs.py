#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 08:50:04 2019

@author: dori
"""

import pandas as pd
import xarray as xr
#import netCDF4 as nc
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from statistic import hist_and_plot
from scipy.optimize import curve_fit
plt.close('all')

def noisePower(h, p0):
  return p0*h**2

filesW = glob('/data/data_hatpro/jue/data/joyrad94/l0/201511/2*/joyrad94_joyce_2015112*.nc')
for fw in filesW:
  print(fw)
  ncfile = xr.open_dataset(fw, drop_variables='velocity')
  Ze = ncfile['Ze']
  Ze = Ze.where(Ze != -999.0)
  Zlin = 10.0**(0.1*Ze)
  t, r = xr.broadcast(Ze.time, Ze.range)
  df = pd.DataFrame()
  df['Hgt'] = r.data.flatten()
  df['Ze'] = Ze.data.flatten()
  spec = ncfile.spec
  spec = spec.where(spec != -999.0)
  speclin = 10.0**(0.1*spec)
  Zg = speclin.sum(dim='velocity', skipna=True)
  N = 10.0*np.log10(Zg - Zlin).data.flatten()
  N[~np.isfinite(N)] = np.nan
  df['N'] = N
  df.dropna(inplace=True, subset=['Ze'])
  df.to_hdf('joyrad94snr.h5', key='stat', mode='a', append=True)

print('done')
df = pd.read_hdf('joyrad94snr.h5', key='stat')
df['SNR'] = df['Ze'] - df['N']
hist_and_plot(df.dropna(subset=['SNR']), 'SNR Joyrad94', yvar='Hgt', xvar='SNR',
              xlabel='SNRW   [dB]', ylabel='Hgt   [m]',
              xlim=[-30, 60], ylim=[0, 12000], lognorm=True,
              savename='SNRw.png', inverty=False, figax=None,
              bins=100, density=False, CFAD=False)

hist_and_plot(df.dropna(subset=['N']), 'Noise Joyrad94', yvar='Hgt', xvar='N',
              xlabel='NoiseW   [dB]', ylabel='Hgt   [m]',
              xlim=[-60, -10], ylim=[0, 12000], lognorm=True,
              savename='NoiseW.png', inverty=False, figax=None,
              bins=100, density=False, CFAD=False)

hist_and_plot(df.dropna(subset=['Ze']), 'Ze Joyrad94', yvar='Hgt', xvar='Ze',
              xlabel='Ze W   [dBZ]', ylabel='Hgt   [m]',
              xlim=[-60, 40], ylim=[0, 12000], lognorm=True,
              savename='ZeW.png', inverty=False, figax=None,
              bins=100, density=False, CFAD=False)

df = df.dropna(subset=['N'])
optW, covW = curve_fit(noisePower, df['Hgt']/1000.0, 10.0**(0.1*df['N']))
print('W pNoise0:', 10.0*np.log10(optW))


#filesK = glob('/data/joyce/data/mmclx/2015/1511/2015112*.gz')
#for fk in filesK:
#  print(fk)
#  ncfile = xr.open_dataset(fk)
#  if (np.min(ncfile.elv)<89.0 or np.max(ncfile.elv)>91):
#    print('skipping scanning data')
#  else:
#    snr = ncfile['SNRg']
#    df = pd.DataFrame()
#    df['SNR'] = 10.0*np.log10(snr.data.flatten())
#    t,r = xr.broadcast(snr.time, snr.range)
#    df['Hgt'] = r.data.flatten()
#    Ze = 10.0*np.log10(ncfile['Zg'].data).flatten()
#    Ze[~np.isfinite(Ze)] = np.nan
#    df['Ze'] = Ze
#    df.dropna(inplace=True, subset=['Ze'])
#    df.to_hdf('joyrad35snr.h5', key='stat', mode='a', append=True)
#
#print('done')
#df = pd.read_hdf('joyrad35snr.h5', key='stat')
#hist_and_plot(df.dropna(subset=['SNR']), 'SNR Joyrad35', yvar='Hgt', xvar='SNR',
#              xlabel='SNRka   [dB]', ylabel='Hgt   [m]',
#              xlim=[-30, 60], ylim=[0, 12000], lognorm=True,
#              savename='SNRka.png', inverty=False, figax=None,
#              bins=100, density=False, CFAD=False)
#
#df['N'] = df['Ze'] - df['SNR']
#hist_and_plot(df.dropna(subset=['N']), 'Noise Joyrad35', yvar='Hgt', xvar='N',
#              xlabel='NoiseKa   [dB]', ylabel='Hgt   [m]',
#              xlim=[-60, -10], ylim=[0, 12000], lognorm=True,
#              savename='NoiseKa.png', inverty=False, figax=None,
#              bins=100, density=False, CFAD=False)
#
#hist_and_plot(df.dropna(subset=['Ze']), 'Ze Joyrad35', yvar='Hgt', xvar='Ze',
#              xlabel='Ze Ka   [dBZ]', ylabel='Hgt   [m]',
#              xlim=[-60, 40], ylim=[0, 12000], lognorm=True,
#              savename='ZeKa.png', inverty=False, figax=None,
#              bins=100, density=False, CFAD=False)
#
#optK, covK = curve_fit(noisePower, df['Hgt']/1000.0, 10.0**(0.1*df['N']))
#print('Ka pNoise0:', 10.0*np.log10(optK))

#df[df['Ze']<-30] = np.nan
#hist_and_plot(df.dropna(subset=['SNR']), 'SNR Joyrad35', yvar='Hgt', xvar='SNR',
#              xlabel='SNRka   [dB]', ylabel='Hgt   [m]',
#              xlim=[-2, 32], ylim=[0, 12000], lognorm=True,
#              savename='SNRkaNoiseCorr.png', inverty=False, figax=None,
#              bins=50, density=False, CFAD=False)
#hist_and_plot(df.dropna(subset=['N']), 'Noise Joyrad35', yvar='Hgt', xvar='N',
#              xlabel='NoiseKa   [dB]', ylabel='Hgt   [m]',
#              xlim=[-60, 0], ylim=[0, 12000], lognorm=True,
#              savename='NoiseKaNoisecorr.png', inverty=False, figax=None,
#              bins=50, density=False, CFAD=False)
#
#hist_and_plot(df.dropna(subset=['Ze']), 'Ze Joyrad35', yvar='Hgt', xvar='Ze',
#              xlabel='Ze Ka   [dBZ]', ylabel='Hgt   [m]',
#              xlim=[-40, 40], ylim=[0, 12000], lognorm=True,
#              savename='ZeKaNoisecorr.png', inverty=False, figax=None,
#              bins=50, density=False, CFAD=False)


#filesX = glob('/data/data_hatpro/jue/tripex/kixpol/*_kixpol.nc')
#i = 0
#chunk = 5
#j = 1
#for fx in filesX:
#  print(fx, i, chunk, str(chunk*j))
#  ncfile = xr.open_dataset(fx)
#  snr = ncfile['SNR']
#  df = pd.DataFrame()
#  df['SNR'] = snr.data.flatten()
#  t,r = xr.broadcast(snr.time, snr.range)
#  df['Hgt'] = r.data.flatten()
#  df['Ze'] = ncfile['Ze'].data.flatten()
#  df.dropna(inplace=True, subset=['Ze'])
#  df.to_hdf('kixpolsnr'+str(chunk*j)+'.h5', key='stat', mode='a', append=True)
#  if i < chunk:
#    i = i + 1
#  else:
#    i = 0
#    j = j + 1

#print('done')
#df = pd.read_hdf('kixpolsnr5.h5', key='stat')
#hist_and_plot(df.dropna(subset=['SNR']), 'SNR KiXPol', yvar='Hgt', xvar='SNR',
#              xlabel='SNRx   [dB]', ylabel='Hgt   [m]',
#              xlim=[-5, 32], ylim=[0, 12000], lognorm=True,
#              savename='SNRx.png', inverty=False, figax=None,
#              bins=50, density=False, CFAD=False)
#
#df['N'] = df['Ze'] - df['SNR']
#hist_and_plot(df.dropna(subset=['N']), 'Noise KiXPol', yvar='Hgt', xvar='N',
#              xlabel='NoiseX   [dB]', ylabel='Hgt   [m]',
#              xlim=[-60, 0], ylim=[0, 12000], lognorm=True,
#              savename='Noisex.png', inverty=False, figax=None,
#              bins=50, density=False, CFAD=False)
#
#hist_and_plot(df.dropna(subset=['Ze']), 'Ze KiXPol', yvar='Hgt', xvar='Ze',
#              xlabel='Ze X   [dBZ]', ylabel='Hgt   [m]',
#              xlim=[-40, 40], ylim=[0, 12000], lognorm=True,
#              savename='ZeX.png', inverty=False, figax=None,
#              bins=50, density=False, CFAD=False)
#
#df[df['Ze']<-30] = np.nan
#hist_and_plot(df.dropna(subset=['SNR']), 'SNR KiXPol', yvar='Hgt', xvar='SNR',
#              xlabel='SNRx   [dB]', ylabel='Hgt   [m]',
#              xlim=[-2, 32], ylim=[0, 12000], lognorm=True,
#              savename='SNRxNoiseCorr.png', inverty=False, figax=None,
#              bins=50, density=False, CFAD=False)
#hist_and_plot(df.dropna(subset=['N']), 'Noise KiXPol', yvar='Hgt', xvar='N',
#              xlabel='NoiseX   [dB]', ylabel='Hgt   [m]',
#              xlim=[-60, 0], ylim=[0, 12000], lognorm=True,
#              savename='NoisexNoisecorr.png', inverty=False, figax=None,
#              bins=50, density=False, CFAD=False)
#
#hist_and_plot(df.dropna(subset=['Ze']), 'Ze KiXPol', yvar='Hgt', xvar='Ze',
#              xlabel='Ze X   [dBZ]', ylabel='Hgt   [m]',
#              xlim=[-40, 40], ylim=[0, 12000], lognorm=True,
#              savename='ZeXNoisecorr.png', inverty=False, figax=None,
#              bins=50, density=False, CFAD=False)
#
#optX, covX = curve_fit(noisePower, df['Hgt']/1000.0, 10.0**(0.1*df['N']))
#print('X pNoise0:', 10.0*np.log10(optX))