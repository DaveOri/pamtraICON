#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 11:06:17 2019

@author: dori
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from statistics import hist_and_plot
from READ import read_prepare
plt.close('all')

addlabel = ''
addlabel = 'run'
#addlabel = 'spin'
hydroset='all_hydro'
data = read_prepare(hydroset, maxhour=46.0, minhour=6.0)
lognorm=True

h,x,y = hist_and_plot(data, 'CFAD  X-band  SNR', 'Hgt', 'N10',
                      'SNR X-band [dBZ]', 'Height   [m]',
                      xlim=[-15, 70], ylim=[0, 10000], inverty=False,
                      savename='CFAD/CFAD_snrX_Hgt' + addlabel +'.png',lognorm=lognorm)

h,x,y = hist_and_plot(data, 'CFAD  Ka-band  SNR', 'Hgt', 'N35',
                      'SNR Ka-band [dB]', 'Height   [m]',
                      xlim=[-15, 70], ylim=[0, 10000], inverty=False,
                      savename='CFAD/CFAD_snrKa_Hgt' + addlabel +'.png', lognorm=lognorm)

h,x,y = hist_and_plot(data, 'CFAD  W-band  SNR', 'Hgt', 'N94',
                      'SNR W-band [dB]', 'Height   [m]',
                      xlim=[-15, 70], ylim=[0, 10000], inverty=False,
                      savename='CFAD/CFAD_snrW_Hgt' + addlabel +'.png',lognorm=lognorm)

#%%############################################################################
data['PN10'] = data['Z10']-data['N10']
data['PN35'] = data['Z35']-data['N35']
data['PN94'] = data['Z94']-data['N94']

def noiseCurve(x, p0):
  return p0 + 20.0*np.log10(x*0.001)

p10, cov10 = curve_fit(noiseCurve, data['Hgt'], data['PN10'])
p35, cov35 = curve_fit(noiseCurve, data['Hgt'], data['PN35'])
p94, cov94 = curve_fit(noiseCurve, data['Hgt'], data['PN94'])
xHgt = np.linspace(100.0,14000.0,1000)

h,x,y = hist_and_plot(data, 'CFAD  X-band  Pnoise', 'Hgt', 'PN10',
                      'Pnoise X-band [dBZ]', 'Height   [m]',
                      xlim=[-70, 0], ylim=[0, 10000], inverty=False,
                      savename='CFAD/CFAD_PnoiseX_Hgt' + addlabel +'.png',lognorm=lognorm)
ax = plt.gca()
ax.plot(noiseCurve(xHgt, *p10), xHgt)
plt.savefig('CFAD/CFAD_PnoiseX_Hgt' + addlabel +'.png')

h,x,y = hist_and_plot(data, 'CFAD  Ka-band  Pnoise', 'Hgt', 'PN35',
                      'Pnoise Ka-band [dBZ]', 'Height   [m]',
                      xlim=[-70, 0], ylim=[0, 10000], inverty=False,
                      savename='CFAD/CFAD_PnoiseKa_Hgt' + addlabel +'.png', lognorm=lognorm)
ax = plt.gca()
ax.plot(noiseCurve(xHgt, *p35), xHgt)
plt.savefig('CFAD/CFAD_PnoiseKa_Hgt' + addlabel +'.png')

h,x,y = hist_and_plot(data, 'CFAD  W-band  Pnoise', 'Hgt', 'PN94',
                      'Pnoise W-band [dBZ]', 'Height   [m]',
                      xlim=[-70, 0], ylim=[0, 10000], inverty=False,
                      savename='CFAD/CFAD_PnoiseW_Hgt' + addlabel +'.png',lognorm=lognorm)
ax = plt.gca()
ax.plot(noiseCurve(xHgt, *p94), xHgt)
plt.savefig('CFAD/CFAD_PnoiseW_Hgt' + addlabel +'.png')

#%%############################################################################

h,x,y = hist_and_plot(data, 'CFAD  X-band  Ze', 'Hgt', 'Z10',
                      'Z X-band [dBZ]', 'Height   [m]',
                      xlim=[-70, 50], ylim=[0, 10000], inverty=False,
                      savename='CFAD/CFAD_X_Hgt' + addlabel +'.png',lognorm=lognorm)
ax = plt.gca()
ax.plot(noiseCurve(xHgt, *p10), xHgt)
plt.savefig('CFAD/CFAD_X_Hgt' + addlabel +'.png')

h,x,y = hist_and_plot(data, 'CFAD  Ka-band  Ze', 'Hgt', 'Z35',
                      'Z Ka-band [dBZ]', 'Height   [m]',
                      xlim=[-70, 50], ylim=[0, 10000], inverty=False,
                      savename='CFAD/CFAD_Ka_Hgt' + addlabel +'.png', lognorm=lognorm)
ax = plt.gca()
ax.plot(noiseCurve(xHgt, *p35), xHgt)
plt.savefig('CFAD/CFAD_Ka_Hgt' + addlabel +'.png')

h,x,y = hist_and_plot(data, 'CFAD  W-band  Ze', 'Hgt', 'Z94',
                      'Z W-band [dBZ]', 'Height   [m]',
                      xlim=[-70, 50], ylim=[0, 10000], inverty=False,
                      savename='CFAD/CFAD_W_Hgt' + addlabel +'.png',lognorm=lognorm)
ax = plt.gca()
ax.plot(noiseCurve(xHgt, *p94), xHgt)
plt.savefig('CFAD/CFAD_W_Hgt' + addlabel +'.png')

#h,x,y = hist_and_plot(data, 'CFAD  DWRxk', 'Hgt', 'DWRxk',
#                      'DWR X-Ka [dB]', 'Height   [m]',
#                      xlim=[-10, 20], ylim=[0, 10000], inverty=False,
#                      savename='CFAD/CFAD_DWRxk_Hgt' + addlabel +'.png', lognorm=lognorm)
#
#h,x,y = hist_and_plot(data, 'CFAD  DWRkw', 'Hgt', 'DWRkw',
#                      'DWR Ka-W [dB]', 'Height   [m]',
#                      xlim=[-10, 20], ylim=[0, 10000], inverty=False,
#                      savename='CFAD/CFAD_DWRkw_Hgt' + addlabel +'.png', lognorm=lognorm)