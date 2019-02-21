# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 11:41:07 2019

@author: dori
"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

dataset = nc.Dataset('/data/obs/site/jue/joyrad35/2019/01/15/' +
                     '20190115_0000.znc')
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

plt.close('all')
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

plt.close('all')
plt.plot(10.0*np.log10(noise_X_lin[0,:]), range_X[0,:])
plt.grid()
plt.xlabel('Noise level X [dB]')
plt.ylabel('range    [m]')
plt.savefig('noise_level_X.png')

plt.figure()
plt.plot(10.0*np.log10(noise_ka_lin[0,:]), range_ka[0,:])
plt.grid()
plt.xlim([-70,-60])
plt.ylim([500,1500])
plt.xlabel('Noise level X [dB]')
plt.ylabel('range    [m]')
plt.savefig('noise_level_X_zoom1k.png')