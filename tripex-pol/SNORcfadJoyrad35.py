#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 17:48:55 2019

@author: dori
"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors

#import sys
#import os
#os.chdir('../comparison/')
#
#from READ import read_prepare
#from READ import read_radar
#from READ import slice_data
#from READ import icon100heights
#from statistics import hist_and_plot
#
#campaign = 'tripex-pol'
#minhour=6.0
#
##hydroset = 'only_ice'
##hydroset = 'only_snow'
##hydroset = 'no_snow'
##hydroset = 'only_graupel_hail'
##hydroset = 'only_liquid'
#hydroset = 'all_hydro'
#
#pamtra = read_prepare(hydroset=campaign + '_' + hydroset + '_', minhour=minhour)
#Hgt = sorted(pamtra['Hgt'].drop_duplicates())
#hist_and_plot(pamtra, 'SNR simulated X', yvar='Hgt', xvar='N10',
#              xlabel='SNR   [dB]', ylabel='Hgt   [m]',
#              xlim=[-40, 60], ylim=[0, 14000], lognorm=True,
#              savename='../tripex-pol/pamtraSNR10.png',
#              inverty=False, figax=None,
#              bins=(100,Hgt))
#hist_and_plot(pamtra, 'SNR simulated Ka', yvar='Hgt', xvar='N35',
#              xlabel='SNR   [dB]', ylabel='Hgt   [m]',
#              xlim=[-40, 60], ylim=[0, 14000], lognorm=True,
#              savename='../tripex-pol/pamtraSNR35.png',
#              inverty=False, figax=None,
#              bins=(100,Hgt))

from glob import glob

def extract(ncfile, varname):
  print(ncfile)
  with nc.Dataset(ncfile) as dataset:
    return dataset.variables[varname][:]

def tiled(ncfile, varname, tile):
  print(ncfile)
  with nc.Dataset(ncfile) as dataset:
    Ntile = len(dataset.variables[tile])
    return np.tile(dataset.variables[varname][:],(Ntile,1))

def dB(value):
  return 10.0*np.log10(value)

def Bd(value):
  return 10.0**(0.1*value)

datafiles = sorted(glob('/data/obs/site/jue/joyrad35/2019/01/1[0-7]/2019011[0-7]_*.znc'))
SNRlist = [extract(df, 'SNRgc') for df in datafiles]
Rlist = [tiled(df, 'range', 'time') for df in datafiles]
SNR = np.concatenate(SNRlist)
R = np.concatenate(Rlist)

H, x, y = np.histogram2d(dB(SNR.flatten()), R.flatten(), bins=100)
xc = (x[:-1] + x[1:])*0.5
yc = (y[:-1] + y[1:])*0.5

plt.figure()
plt.pcolormesh(xc, yc, H, norm=colors.LogNorm(vmin=1, vmax=np.nanmax(H)),
               cmap='jet')
plt.xlabel('Joyrad35 SNRg')
plt.ylabel('Height  [m]')
plt.title('Joyrad35')
plt.savefig('SNRjrad35.png')

datafiles = sorted(glob('/data/obs/site/jue/joyrad10/2019/01/1[0-7]/2019011[0-7]_*.znc'))
SNRlist = [extract(df, 'SNRgc') for df in datafiles]
Rlist = [tiled(df, 'range', 'time') for df in datafiles]
SNR10 = np.concatenate(SNRlist)
R10 = np.concatenate(Rlist)

H10, x10, y10 = np.histogram2d(dB(SNR10.flatten()), R10.flatten(), bins=100)
xc10 = (x10[:-1] + x10[1:])*0.5
yc10 = (y10[:-1] + y10[1:])*0.5

plt.figure()
plt.pcolormesh(xc10, yc10, H10, norm=colors.LogNorm(vmin=1, vmax=np.nanmax(H)),
               cmap='jet')
plt.xlabel('Joyrad10 SNRg')
plt.ylabel('Height  [m]')
plt.title('Joyrad10')
plt.savefig('SNRjrad10.png')