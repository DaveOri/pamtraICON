#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 14:19:33 2019

@author: dori
"""

#import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
from statistics import hist_and_plot
from READ import read_prepare
plt.close('all')

addlabel = ''
addlabel = 'run'
#addlabel = 'spin'
campaign = 'tripex'
hydroset = campaign + '_all_hydro_'
data = read_prepare(hydroset, maxhour=6.0, minhour=0.0)

lognorm=True
h,x,y = hist_and_plot(data, 'CFAD  X-band  Ze', 'T', 'Z10',
                      'Z X-band [dBZ]', 'Temperature [deg C]',
                      xlim=[-70, 50], ylim=[-60, 20], inverty=True,
                      savename='CFAD/CFAD_X_temp' + campaign  + addlabel +'.png',
                      lognorm=lognorm)

h,x,y = hist_and_plot(data, 'CFAD  Ka-band  Ze', 'T', 'Z35',
                      'Z Ka-band [dBZ]', 'Temperature [deg C]',
                      xlim=[-70, 50], ylim=[-60, 20], inverty=True,
                      savename='CFAD/CFAD_Ka_temp' + campaign  + addlabel +'.png',
                      lognorm=lognorm)

h,x,y = hist_and_plot(data, 'CFAD  W-band  Ze', 'T', 'Z94',
                      'Z W-band [dBZ]', 'Temperature [deg C]',
                      xlim=[-70, 50], ylim=[-60, 20], inverty=True,
                      savename='CFAD/CFAD_W_temp' + campaign  + addlabel +'.png',
                      lognorm=lognorm)

h,x,y = hist_and_plot(data, 'CFAD  DWRxk', 'T', 'DWRxk',
                      'DWR X-Ka [dB]', 'Temperature [deg C]',
                      xlim=[-10, 20], ylim=[-60, 20], inverty=True,
                      savename='CFAD/CFAD_DWRxk_temp' + campaign  + addlabel +'.png',
                      lognorm=lognorm)

h,x,y = hist_and_plot(data, 'CFAD  DWRkw', 'T', 'DWRkw',
                      'DWR Ka-W [dB]', 'Temperature [deg C]',
                      xlim=[-10, 20], ylim=[-60, 20], inverty=True,
                      savename='CFAD/CFAD_DWRkw_temp' + campaign  + addlabel +'.png',
                      lognorm=lognorm)

#%%############################################################################
lognorm=True
h,x,y = hist_and_plot(data, 'CFAD  X-band  MDV [m/s]', 'T', 'V10',
                      'MDV X-band [m/s]', 'T   [deg C]',
                      xlim=[-8, 2], ylim=[-60, 20], inverty=True,
                      savename='CFAD/CFAD_velX_temp' + campaign  + addlabel +'.png',
                      lognorm=lognorm)

h,x,y = hist_and_plot(data, 'CFAD  Ka-band  MDV [m/s]', 'T', 'V35',
                      'MDV Ka-band [m/s]', 'T   [deg C]',
                      xlim=[-8, 2], ylim=[-60, 20], inverty=True,
                      savename='CFAD/CFAD_velKa_temp' + campaign  + addlabel +'.png',
                      lognorm=lognorm)

h,x,y = hist_and_plot(data, 'CFAD  W-band  MDV [m/s]', 'T', 'V94',
                      'MDV W-band [m/s]', 'T   [deg C]',
                      xlim=[-8, 2], ylim=[-60, 20], inverty=True,
                      savename='CFAD/CFAD_velW_temp' + campaign  + addlabel +'.png',
                      lognorm=lognorm)
#%%############################################################################
lognorm=True
h,x,y = hist_and_plot(data, 'CFAD  X-band  SW [m/s]', 'T', 'W10',
                      'Spectral Width X-band [m/s]', 'T   [deg C]',
                      xlim=[0, 3], ylim=[-60, 20], inverty=True,
                      savename='CFAD/CFAD_widX_temp' + campaign  + addlabel +'.png',
                      lognorm=lognorm)

h,x,y = hist_and_plot(data, 'CFAD  Ka-band  SW [m/s]', 'T', 'W35',
                      'Spectral Width Ka-band [m/s]', 'T   [deg C]',
                      xlim=[0, 3], ylim=[-60, 20], inverty=True,
                      savename='CFAD/CFAD_widKa_temp' + campaign  + addlabel +'.png',
                      lognorm=lognorm)

h,x,y = hist_and_plot(data, 'CFAD  W-band  SW [m/s]', 'T', 'W94',
                      'Spectral Width W-band [m/s]', 'T   [deg C]',
                      xlim=[0, 3], ylim=[-60, 20], inverty=True,
                      savename='CFAD/CFAD_widW_temp' + campaign  + addlabel +'.png',
                      lognorm=lognorm)