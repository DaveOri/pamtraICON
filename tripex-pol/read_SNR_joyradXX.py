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
dataset = nc.Dataset('/data/obs/site/jue/joyrad35/2019/01/13/' +
                     '20190113_0800.znc')
ti, ri = 1430, 23
range_ka = dataset.variables['range'][:]
HSDco = dataset.variables['HSDco'][ti, :]
SNRCorFaCo = dataset.variables['SNRCorFaCo'][ti, :]
radar_const = dataset.variables['RadarConst'][ti]
noise_power_co = dataset.variables['npw1'][ti]
DSP2lin = radar_const*SNRCorFaCo*(range_ka/5000.)**2./noise_power_co
noise_ka_lin = DSP2lin*HSDco

SpectraDSP = dataset.variables['SPCco'][ti, :, :]
Nr, Nd = SpectraDSP.shape
SpectraLin = (SpectraDSP.T*DSP2lin).T
SpectradB = 10.0*np.log10(SpectraLin)
dopplerV = dataset.variables['doppler'][:]


plt.close('all')
plt.figure()
plt.plot(dopplerV, SpectradB[ri, :],
         label='H= {0:.0f}'.format(range_ka[ri]))
plt.plot(dopplerV, 10.0*np.log10(noise_ka_lin[ri])*np.ones(dopplerV.shape),
         label='mean noise lvl')
plt.legend()
plt.grid()

plt.figure()
plt.plot(10.0*np.log10(noise_ka_lin[:]), range_ka[:])
plt.grid()
plt.xlabel('Noise level Ka [dB]')
plt.ylabel('range    [m]')
plt.savefig('noise_level_Ka.png')

plt.figure()
plt.plot(10.0*np.log10(noise_ka_lin[:]), range_ka[:])
plt.grid()
plt.xlim([-70,-60])
plt.ylim([500,1500])
plt.xlabel('Noise level Ka [dB]')
plt.ylabel('range    [m]')
plt.savefig('noise_level_Ka_zoom1k.png')

plt.figure()
plt.plot(10.0*np.log10(noise_ka_lin[:]), range_ka[:])
plt.grid()
plt.xlim([-50,-40])
plt.ylim([8500,9500])
plt.xlabel('Noise level Ka [dB]')
plt.ylabel('range    [m]')
plt.savefig('noise_level_Ka_zoom9k.png')


#%% DO the same for X band Joyrad10
dataset = nc.Dataset('/data/obs/site/jue/joyrad10/2019/01/15/' +
                     '20190115_120004.znc')
dataset = nc.Dataset('/data/obs/site/jue/joyrad10/2019/01/13/' +
                     '20190113_080004.znc')
tt, rr = 1430, 25
range_x = dataset.variables['range'][:]
HSDco = dataset.variables['HSDco'][tt, :]
SNRCorFaCo = dataset.variables['SNRCorFaCo'][tt, :]
radar_const = dataset.variables['RadarConst'][tt]
noise_power_co = dataset.variables['npw1'][tt]
DSP2lin = radar_const*SNRCorFaCo*(range_x/5000.)**2./noise_power_co
noise_x_lin = DSP2lin*HSDco

SpectraDSP = dataset.variables['SPCco'][tt, :, :]
Nr, Nd = SpectraDSP.shape
SpectraLin = SpectraDSP*DSP2lin[:, np.newaxis]#(SpectraDSP.T*DSP2lin).T
SpectradB = 10.0*np.log10(SpectraLin)
dopplerV = dataset.variables['doppler'][:]

plt.close('all')
plt.figure()
plt.plot(dopplerV, SpectradB[rr, :],
         label='H= {0:.0f}'.format(range_x[rr]))
plt.plot(dopplerV, 10.0*np.log10(noise_x_lin[rr])*np.ones(dopplerV.shape),
         label='mean noise lvl')
plt.legend()
plt.grid()

plt.figure()
plt.plot(10.0*np.log10(noise_x_lin[:]), range_x[:])
plt.grid()
plt.xlabel('Noise level X [dB]')
plt.ylabel('range    [m]')
plt.savefig('noise_level_X.png')

plt.figure()
plt.plot(10.0*np.log10(noise_x_lin[:]), range_x[:])
plt.grid()
plt.xlim([-70,-60])
plt.ylim([500,1500])
plt.xlabel('Noise level X [dB]')
plt.ylabel('range    [m]')
plt.savefig('noise_level_X_zoom1k.png')

plt.figure()
plt.plot(10.0*np.log10(noise_x_lin[:]), range_x[:])
plt.grid()
plt.xlim([-50,-40])
plt.ylim([8500,9500])
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