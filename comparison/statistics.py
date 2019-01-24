#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 14:13:11 2019

@author: dori
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.close('all')

data = pd.read_hdf('data_pamtra_icon.h5',key='stat')
data[data==-9999] = np.nan
dataZ = data.dropna()

dataZ['DWRxk'] = dataZ.loc[:,'Z10'] - dataZ.loc[:,'Z35']
dataZ['DWRkw'] = dataZ.loc[:,'Z35'] - dataZ.loc[:,'Z94']

def hist_and_plot(data, title):
  hst, xedge, yedge = np.histogram2d(data['DWRxk'], data['DWRkw'], bins=100)
  xcenter = (xedge[:-1] + xedge[1:])*0.5
  ycenter = (yedge[:-1] + yedge[1:])*0.5
  hst[hst < 1] = np.nan
  fig, ax = plt.subplots(1,1)
  mesh = ax.pcolormesh(xcenter, ycenter, hst[:], cmap='jet',
                       norm=colors.LogNorm(vmin=1, vmax=np.nanmax(hst)))
  ax.set_xlim([-5,20])
  ax.set_ylim([-5,20])
  plt.colorbar(mesh, ax=ax, extend='max', label='counts')
  ax.set_title(title)
  ax.set_xlabel('DWR Ka W   [dB]')
  ax.set_ylabel('DWR X Ka   [dB]')
  ax.grid()
  fig.savefig('triple_frequency/' + '_'.join(title.split()) + '.png')
  return hst, xcenter, ycenter

title = 'contour_log'
fig, ax = plt.subplots(1,1)

T0 = 273.15

def add_one_contour(ax, data, Tminmax, color, levels=1):
  if ax is None:
    fig, ax = plt.subplots(1,1)
  Tmin, Tmax = Tminmax
  if (Tmin is None) and (Tmax is None):
    Tstr = 'whole dataset'
    h,x,y = hist_and_plot(data[:], Tstr)
  elif Tmin is None:
    Tstr = 'T < '+str(Tmax)
    h,x,y = hist_and_plot(data[(data['T'] < Tmax+T0)], Tstr)
  elif Tmax is None:
    Tstr = 'T > '+str(Tmin)
    h,x,y = hist_and_plot(data[(data['T'] > Tmin+T0)], Tstr)
  else:
    Tstr = str(Tmin)+' < T < '+str(Tmax)
    h,x,y = hist_and_plot(data[(data['T'] <= Tmax+T0)*(data['T'] > Tmin+T0)],
                          Tstr)
  CS = ax.contour(x, y, np.log10(h[:]), 1, colors=color)
  CS.collections[0].set_label(Tstr)
  return ax

add_one_contour(None, dataZ, [None, None], 'k')
add_one_contour(None, dataZ, [None, 0], 'k')
add_one_contour(ax, dataZ, [0, None], 'C1')
add_one_contour(ax, dataZ, [-5, 0], 'C2')
add_one_contour(ax, dataZ, [-10, -5], 'C3')
add_one_contour(ax, dataZ, [-15, -10], 'C4')
add_one_contour(ax, dataZ, [-20, -15], 'C5')
add_one_contour(ax, dataZ, [-25, -20], 'C6')
add_one_contour(ax, dataZ, [-30, -25], 'C7')


ax.set_xlim([-5,15])
ax.set_ylim([-5,15])
ax.set_title(title)
ax.set_xlabel('DWR Ka W   [dB]')
ax.set_ylabel('DWR X Ka   [dB]')
ax.grid()
ax.legend()
fig.savefig('triple_frequency/' + '_'.join(title.split()) + '.png')