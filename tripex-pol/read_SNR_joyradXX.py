# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 11:41:07 2019

@author: dori
"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

# Calculate noise level tripex-pol Ka-band Joyrad35
dataset = nc.Dataset('/data/obs/site/jue/joyrad35/2019/01/15/' +
                     '20190115_1200.znc')
range_ka = dataset.variables['range'][:]
HSDco = dataset.variables['HSDco'][:]
Nt, Nr = HSDco.shape
range_ka = np.tile(range_ka,(Nt,1))
SNRCorFaCo = dataset.variables['SNRCorFaCo'][:]
radar_const = dataset.variables['RadarConst'][:]
radar_const = np.tile(radar_const[np.newaxis].T,(1,Nr))
noise_power_co = dataset.variables['npw1'][:]
noise_power_co = np.tile(noise_power_co[np.newaxis].T,(1,Nr))
noise_ka_lin = radar_const*SNRCorFaCo*HSDco/noise_power_co*(range_ka/5000.)**2.

SpectraDSP = dataset.variables['SPCco'][:]
Nt, Nr, Nd = SpectraDSP.shape
radar_const3d = np.tile(radar_const[...,np.newaxis],(1,1,Nd))
range_ka3d = np.tile(range_ka[...,np.newaxis],(1,1,Nd))
SpectradB = 10.0*np.log10(radar_const3d*SpectraDSP*(range_ka3d/5000.)**2)
SpectradB = SpectradB - 10.0*np.log10(dataset.variables['nfft'][0])
dopplerV = dataset.variables['doppler'][:]

ti, ri = 1430, 23
plt.close('all')
plt.figure()
plt.plot(dopplerV, SpectradB[ti, ri, :],
         label='H= {0:.0f}'.format(range_ka[ti, ri]))
plt.plot(dopplerV, 10.0*np.log10(noise_ka_lin[ti, ri])*np.ones(dopplerV.shape),
         label='mean noise lvl')
plt.legend()
plt.grid()

plt.figure()
plt.plot(10.0*np.log10(noise_ka_lin[0,:]), range_ka[0,:])
plt.grid()
plt.xlabel('Noise level Ka [dB]')
plt.ylabel('range    [m]')
plt.savefig('noise_level_Ka.png')

plt.figure()
plt.plot(10.0*np.log10(noise_ka_lin[0,:]), range_ka[0,:])
plt.grid()
plt.xlim([-70,-60])
plt.ylim([500,1500])
plt.xlabel('Noise level Ka [dB]')
plt.ylabel('range    [m]')
plt.savefig('noise_level_Ka_zoom1k.png')

plt.figure()
plt.plot(10.0*np.log10(noise_ka_lin[0,:]), range_ka[0,:])
plt.grid()
plt.xlim([-50,-40])
plt.ylim([8500,9500])
plt.xlabel('Noise level Ka [dB]')
plt.ylabel('range    [m]')
plt.savefig('noise_level_Ka_zoom9k.png')


# DO the same for X band Joyrad10
dataset = nc.Dataset('/data/obs/site/jue/joyrad10/2019/01/15/' +
                     '20190115_000004.znc')
range_X = dataset.variables['range'][:]
HSDco = dataset.variables['HSDco'][:]
Nt, Nr = HSDco.shape
range_X = np.tile(range_X,(Nt,1))
SNRCorFaCo = dataset.variables['SNRCorFaCo'][:]
radar_const = dataset.variables['RadarConst'][:]
radar_const = np.tile(radar_const[np.newaxis].T,(1,Nr))
noise_power_co = dataset.variables['npw1'][:]
noise_power_co = np.tile(noise_power_co[np.newaxis].T,(1,Nr))

noise_X_lin = radar_const*SNRCorFaCo*HSDco/noise_power_co*(range_X/5000.)**2.

#tt, rr = 2778, 102
#SpectraDSP = dataset.variables['SPCco'][tt,rr,:]
#Nt, Nr, Nd = 1,1,1#SpectraDSP.shape
##radar_const3d = np.tile(radar_const[...,np.newaxis],(1,1,Nd))
#radar_const3d = radar_const[tt,rr]
##range_X3d = np.tile(range_X[...,np.newaxis],(1,1,Nd))
#range_X3d = range_X[tt,rr]
#SpectradB = 10.0*np.log10(radar_const3d*SpectraDSP*(range_X3d/5000.)**2)
#dopplerV = dataset.variables['doppler'][:]
#plt.close('all')
#plt.figure()
##plt.plot(dopplerV, SpectradB[tt,rr,:])
#plt.plot(dopplerV, SpectradB)

plt.figure()
plt.plot(10.0*np.log10(noise_X_lin[0,:]), range_X[0,:])
plt.grid()
plt.xlabel('Noise level X [dB]')
plt.ylabel('range    [m]')
plt.savefig('noise_level_X.png')

plt.figure()
plt.plot(10.0*np.log10(noise_X_lin[0,:]), range_X[0,:])
plt.grid()
plt.xlim([-55,-45])
plt.ylim([500,1500])
plt.xlabel('Noise level X [dB]')
plt.ylabel('range    [m]')
plt.savefig('noise_level_X_zoom1k.png')

plt.figure()
plt.plot(10.0*np.log10(noise_X_lin[0,:]), range_X[0,:])
plt.grid()
plt.xlim([-30,-25])
plt.ylim([8000,10000])
plt.xlabel('Noise level X [dB]')
plt.ylabel('range    [m]')
plt.savefig('noise_level_X_zoom9k.png')


#%% Calculate noise level TRIPEx Ka-band Joyrad35
#import gzip
#filename = '/data/joyce/data/mmclx/2015/1511/20151124_1000.mmclx.00.gz'
#with gzip.open(filename) as f:
#  mmclx = f.read()
#ncfile = '/run/user/30376/' + filename.split('/')[-1][:-3]
#with open(ncfile, 'wb') as f:
#  f.write(mmclx)
#dataset = nc.Dataset(ncfile)
#range_ka = dataset.variables['range'][:]
#HSDco = dataset.variables['HSDco'][:]
#Nt, Nr = HSDco.shape
#range_ka = np.tile(range_ka,(Nt,1))
#SNRCorFaCo = dataset.variables['SNRCorFaCo'][:]
#radar_const = dataset.variables['RadarConst'][:]
#radar_const = np.tile(radar_const[np.newaxis].T,(1,Nr))
#noise_power_co = dataset.variables['npw1'][:]
#noise_power_co = np.tile(noise_power_co[np.newaxis].T,(1,Nr))
#
#noise_ka_lin = radar_const*SNRCorFaCo*HSDco/noise_power_co*(range_ka/5000.)**2.
#
#plt.close('all')
#plt.plot(10.0*np.log10(noise_ka_lin[0,:]), range_ka[0,:])
#plt.grid()
#plt.xlabel('Noise level Ka [dB]')
#plt.ylabel('range    [m]')
#plt.savefig('noise_level_Ka_tripex.png')
#
#plt.figure()
#plt.plot(10.0*np.log10(noise_ka_lin[0,:]), range_ka[0,:])
#plt.grid()
#plt.xlim([-70,-60])
#plt.ylim([500,1500])
#plt.xlabel('Noise level Ka [dB]')
#plt.ylabel('range    [m]')
#plt.savefig('noise_level_Ka_zoom1k_tripex.png')