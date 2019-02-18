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
#addlabel = 'spin'
hydroset='only_ice'
data = read_prepare(hydroset, maxhour=40.0, minhour=6.0)

lognorm=True
xlim = [-60, 40]
ylim = [0, 1]
h,x,y = hist_and_plot(data, hydroset + ' SW Z X-band', 'W10', 'Z10',
                      'Z X-band [dBZ]', 'SW  X-band  [m/s]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='SW_Z/SW_Z_X' + addlabel + '_' + hydroset + '.png', lognorm=lognorm)

h,x,y = hist_and_plot(data, hydroset + ' SW Z Ka-band', 'W35', 'Z35',
                      'Z Ka-band [dBZ]', 'SW  Ka-band  [m/s]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='SW_Z/SW_Z_Ka' + addlabel + '_' + hydroset + '.png', lognorm=lognorm)

h,x,y = hist_and_plot(data, hydroset + ' SW Z W-band', 'W94', 'Z94',
                      'Z W-band [dBZ]', 'SW  W-band  [m/s]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='SW_Z/SW_Z_W' + addlabel + '_' + hydroset + '.png', lognorm=lognorm)

hydroset='only_snow'
data = read_prepare(hydroset, maxhour=40.0, minhour=6.0)

lognorm=True
xlim = [-60, 40]
ylim = [0, 1]
h,x,y = hist_and_plot(data, hydroset + ' SW Z X-band', 'W10', 'Z10',
                      'Z X-band [dBZ]', 'SW  X-band  [m/s]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='SW_Z/SW_Z_X' + addlabel + '_' + hydroset + '.png', lognorm=lognorm)

h,x,y = hist_and_plot(data, hydroset + ' SW Z Ka-band', 'W35', 'Z35',
                      'Z Ka-band [dBZ]', 'SW  Ka-band  [m/s]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='SW_Z/SW_Z_Ka' + addlabel + '_' + hydroset + '.png', lognorm=lognorm)

h,x,y = hist_and_plot(data, hydroset + ' SW Z W-band', 'W94', 'Z94',
                      'Z W-band [dBZ]', 'SW  W-band  [m/s]',
                      xlim=xlim, ylim=ylim, inverty=False,
                      savename='SW_Z/SW_Z_W' + addlabel + '_' + hydroset + '.png', lognorm=lognorm)