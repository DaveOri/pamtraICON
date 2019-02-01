#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 11:21:10 2019

@author: dori
"""

#import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
from statistics import hist_and_plot, make6panel_temperature
from READ import read_prepare
plt.close('all')

addlabel = ''
addlabel = 'run'
addlabel = 'spin'
hydroset='all_hydro'
data = read_prepare(hydroset, maxhour=6.0, minhour=0.0)

lognorm=True
xlim = [-10, 8]
ylim = [-5, 20]
h,x,y = hist_and_plot(data, 'DWRxk MDV X-band', 'DWRxk', 'V10',
                      'MDV X-band [m/s]', 'DWRxk [dB]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='DWRvel/DWRxk_MDV_X' + addlabel +'.png', lognorm=lognorm)

h,x,y = hist_and_plot(data, 'DWRxk MDV Ka-band', 'DWRxk', 'V35',
                      'MDV Ka-band [m/s]', 'DWRxk [dB]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='DWRvel/DWRxk_MDV_Ka' + addlabel +'.png', lognorm=lognorm)

h,x,y = hist_and_plot(data, 'DWRxk MDV W-band', 'DWRxk', 'V94',
                      'MDV W-band [m/s]', 'DWRxk [dB]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='DWRvel/DWRxk_MDV_W' + addlabel +'.png', lognorm=lognorm)

h,x,y = hist_and_plot(data, 'DWRkw MDV X-band', 'DWRkw', 'V10',
                      'MDV X-band [m/s]', 'DWRkw [dB]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='DWRvel/DWRkw_MDV_X' + addlabel +'.png', lognorm=lognorm)

h,x,y = hist_and_plot(data, 'DWRkw MDV Ka-band', 'DWRkw', 'V35',
                      'MDV Ka-band [m/s]', 'DWRkw [dB]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='DWRvel/DWRkw_MDV_Ka' + addlabel +'.png', lognorm=lognorm)

h,x,y = hist_and_plot(data, 'DWRkw MDV W-band', 'DWRkw', 'V94',
                      'MDV W-band [m/s]', 'DWRkw [dB]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='DWRvel/DWRkw_MDV_W' + addlabel +'.png', lognorm=lognorm)

#%%############################################################################

plt.close('all')
lognorm=True
xlim = [0, 3]
ylim = [-5, 20]

h,x,y = hist_and_plot(data, 'DWRxk Spec Width X-band', 'DWRxk', 'W10',
                      'Spec Width X-band [m/s]', 'DWRxk [dB]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='DWRvel/DWRxk_SW_X' + addlabel +'.png', lognorm=lognorm)

h,x,y = hist_and_plot(data, 'DWRxk Spec Width Ka-band', 'DWRxk', 'W35',
                      'Spec Width Ka-band [m/s]', 'DWRxk [dB]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='DWRvel/DWRxk_SW_Ka' + addlabel +'.png', lognorm=lognorm)

h,x,y = hist_and_plot(data, 'DWRxk Spec Width W-band', 'DWRxk', 'W94',
                      'Spec Width W-band [m/s]', 'DWRxk [dB]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='DWRvel/DWRxk_SW_W' + addlabel +'.png', lognorm=lognorm)

h,x,y = hist_and_plot(data, 'DWRkw Spec Width X-band', 'DWRkw', 'W10',
                      'Spec Width X-band [m/s]', 'DWRkw [dB]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='DWRvel/DWRkw_SW_X' + addlabel +'.png', lognorm=lognorm)

h,x,y = hist_and_plot(data, 'DWRkw Spec Width Ka-band', 'DWRkw', 'W35',
                      'Spec Width Ka-band [m/s]', 'DWRkw [dB]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='DWRvel/DWRkw_SW_Ka' + addlabel +'.png', lognorm=lognorm)

h,x,y = hist_and_plot(data, 'DWRkw Spec Width W-band', 'DWRkw', 'W94',
                      'Spec Width W-band [m/s]', 'DWRkw [dB]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='DWRvel/DWRkw_SW_W' + addlabel +'.png', lognorm=lognorm)

#%%############################################################################
xlim = [-10, 2]
ylim = [-5, 20]
fig, axs = plt.subplots(3,2, sharex=True, figsize=[12.8,9.6])
h,x,y = hist_and_plot(data[data['T'] > 0],
                      'T > 0', 'DWRxk', 'V10',
                      None, 'DWRxk [dB]', figax=[fig, axs[0,0]],
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename=None, lognorm=lognorm)

h,x,y = hist_and_plot(data[(data['T'] > -5)*(data['T'] < -0)],
                      '-5 < T < 0', 'DWRxk', 'V10',
                      None, 'DWRxk [dB]', figax=[fig, axs[0,1]],
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename=None, lognorm=lognorm)

h,x,y = hist_and_plot(data[(data['T'] > -10)*(data['T'] < -5)],
                      '-10 < T < -5', 'DWRxk', 'V10',
                      None, 'DWRxk [dB]', figax=[fig, axs[1,0]],
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename=None, lognorm=lognorm)

h,x,y = hist_and_plot(data[(data['T'] > -15)*(data['T'] < -10)],
                      '-15 < T < -10', 'DWRxk', 'V10',
                      None, 'DWRxk [dB]', figax=[fig, axs[1,1]],
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename=None, lognorm=lognorm)

h,x,y = hist_and_plot(data[(data['T'] > -20)*(data['T'] < -15)],
                      '-20 < T < -15', 'DWRxk', 'V10',
                      'MDV X-band [m/s]', 'DWRxk [dB]', figax=[fig, axs[2,0]],
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename=None, lognorm=lognorm)

h,x,y = hist_and_plot(data[(data['T'] > -25)*(data['T'] < -20)],
                      '-25 < T < -20', 'DWRxk', 'V10',
                      'MDV X-band [m/s]', 'DWRxk [dB]', figax=[fig, axs[2,1]],
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename=None, lognorm=lognorm)
fig.savefig('DWRvel/DWRxk_MDV_X_temp' + addlabel +'.png')

#%%############################################################################
xlim = [0, 3]
ylim = [-5, 20]

make6panel_temperature(data, 'DWRxk', 'W10', xlim, ylim, lognorm,
                       'Spec Width X-band [m/s]', 'DWRxk [dB]',
                       'DWRvel/DWRxk_SW_X_temp' + addlabel +'.png')

make6panel_temperature(data, 'DWRxk', 'W35', xlim, ylim, lognorm,
                       'Spec Width Ka-band [m/s]', 'DWRxk [dB]',
                       'DWRvel/DWRxk_SW_Ka_temp' + addlabel +'.png')

make6panel_temperature(data, 'DWRxk', 'W94', xlim, ylim, lognorm,
                       'Spec Width W-band [m/s]', 'DWRxk [dB]',
                       'DWRvel/DWRxk_SW_W_temp' + addlabel +'.png')

make6panel_temperature(data, 'DWRkw', 'W10', xlim, ylim, lognorm,
                       'Spec Width X-band [m/s]', 'DWRkw [dB]',
                       'DWRvel/DWRkw_SW_X_temp' + addlabel +'.png')

make6panel_temperature(data, 'DWRkw', 'W35', xlim, ylim, lognorm,
                       'Spec Width Ka-band [m/s]', 'DWRkw [dB]',
                       'DWRvel/DWRkw_SW_Ka_temp' + addlabel +'.png')

make6panel_temperature(data, 'DWRkw', 'W94', xlim, ylim, lognorm,
                       'Spec Width W-band [m/s]', 'DWRkw [dB]',
                       'DWRvel/DWRkw_SW_W_temp' + addlabel +'.png')

xlim = [-10, 2]
ylim = [-5, 20]

make6panel_temperature(data, 'DWRxk', 'V10', xlim, ylim, lognorm,
                       'MDV X-band [m/s]', 'DWRxk [dB]',
                       'DWRvel/DWRxk_MDV_X_temp' + addlabel +'.png')

make6panel_temperature(data, 'DWRxk', 'V35', xlim, ylim, lognorm,
                       'MDV Ka-band [m/s]', 'DWRxk [dB]',
                       'DWRvel/DWRxk_MDV_Ka_temp' + addlabel +'.png')

make6panel_temperature(data, 'DWRxk', 'V94', xlim, ylim, lognorm,
                       'MDV W-band [m/s]', 'DWRxk [dB]',
                       'DWRvel/DWRxk_MDV_W_temp' + addlabel +'.png')

make6panel_temperature(data, 'DWRkw', 'V10', xlim, ylim, lognorm,
                       'MDV X-band [m/s]', 'DWRkw [dB]',
                       'DWRvel/DWRkw_MDV_X_temp' + addlabel +'.png')

make6panel_temperature(data, 'DWRkw', 'V35', xlim, ylim, lognorm,
                       'MDV Ka-band [m/s]', 'DWRkw [dB]',
                       'DWRvel/DWRkw_MDV_Ka_temp' + addlabel +'.png')

make6panel_temperature(data, 'DWRkw', 'V94', xlim, ylim, lognorm,
                       'MDV W-band [m/s]', 'DWRkw [dB]',
                       'DWRvel/DWRkw_MDV_W_temp' + addlabel +'.png')
