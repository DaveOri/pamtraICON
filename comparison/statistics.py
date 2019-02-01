#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 14:13:11 2019

@author: dori
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from READ import read_prepare
plt.close('all')

addlabel = ''
addlabel = 'run'
#addlabel = 'spin'
hydroset='all_hydro'
data = read_prepare(hydroset, maxhour=46.0, minhour=9.0)

def hist_and_plot(data, title, yvar='DWRxk', xvar='DWRkw',
                  xlabel='DWR Ka W   [dB]', ylabel='DWR X Ka   [dB]',
                  xlim=[-5, 20], ylim=[-5, 20], lognorm=True,
                  savename='auto3f', inverty=False, figax=None):
  hst, xedge, yedge = np.histogram2d(data[xvar], data[yvar], bins=100)
  hst = hst.T
  xcenter = (xedge[:-1] + xedge[1:])*0.5
  ycenter = (yedge[:-1] + yedge[1:])*0.5
  hst[hst < 1] = np.nan
  if figax is None:
    fig, ax = plt.subplots(1,1)
  else:
    fig, ax = figax
  if lognorm:
    mesh = ax.pcolormesh(xcenter, ycenter, hst[:], cmap='jet',
                         norm=colors.LogNorm(vmin=1, vmax=np.nanmax(hst)))
  else:
    mesh = ax.pcolormesh(xcenter, ycenter, hst[:], cmap='jet')
  ax.set_xlim(xlim)
  ax.set_ylim(ylim)
  plt.colorbar(mesh, ax=ax, extend='max', label='counts')
  ax.set_title(title)
  if xlabel is not None:
    ax.set_xlabel(xlabel)
  if ylabel is not None:
    ax.set_ylabel(ylabel)
  ax.grid()
  if inverty:
    ax.invert_yaxis()
  if savename is 'auto3f':
    fig.savefig('triple_frequency/' + '_'.join(title.split()) + addlabel + '.png')
  elif savename is None:
    pass
  else:
    fig.savefig(savename)
  return hst, xcenter, ycenter

def add_one_contour(ax, data, Tminmax, color, levels=1):
  if ax is None:
    fig, ax = plt.subplots(1,1)
  Tmin, Tmax = Tminmax
  if (Tmin is None) and (Tmax is None):
    Tstr = 'whole dataset'
    h,x,y = hist_and_plot(data[:], Tstr)
  elif Tmin is None:
    Tstr = 'T < '+str(Tmax)
    h,x,y = hist_and_plot(data[(data['T'] < Tmax)], Tstr)
  elif Tmax is None:
    Tstr = 'T > '+str(Tmin)
    h,x,y = hist_and_plot(data[(data['T'] > Tmin)], Tstr)
  else:
    Tstr = str(Tmin)+' < T < '+str(Tmax)
    h,x,y = hist_and_plot(data[(data['T'] <= Tmax)*(data['T'] > Tmin)],
                          Tstr)
  CS = ax.contour(x, y, np.log10(h[:]), 1, colors=color)
  CS.collections[0].set_label(Tstr)
  return ax

def make6panel_temperature(data, yvar, xvar, xlim, ylim, lognorm, xlab, ylab, savename):
  fig, axs = plt.subplots(3,2, sharex=True, figsize=[12.8,9.6])
  h,x,y = hist_and_plot(data[data['T'] > 0],
                        'T > 0', yvar, xvar,
                        None, ylab, figax=[fig, axs[0,0]],
                        xlim=xlim, ylim=ylim, inverty=False,
                        savename=None, lognorm=lognorm)

  h,x,y = hist_and_plot(data[(data['T'] > -5)*(data['T'] < -0)],
                        '-5 < T < 0', yvar, xvar,
                        None, ylab, figax=[fig, axs[0,1]],
                        xlim=xlim, ylim=ylim, inverty=False,
                        savename=None, lognorm=lognorm)
  
  h,x,y = hist_and_plot(data[(data['T'] > -10)*(data['T'] < -5)],
                        '-10 < T < -5', yvar, xvar,
                        None, ylab, figax=[fig, axs[1,0]],
                        xlim=xlim, ylim=ylim, inverty=False,
                        savename=None, lognorm=lognorm)
  
  h,x,y = hist_and_plot(data[(data['T'] > -15)*(data['T'] < -10)],
                        '-15 < T < -10', yvar, xvar,
                        None, ylab, figax=[fig, axs[1,1]],
                        xlim=xlim, ylim=ylim, inverty=False,
                        savename=None, lognorm=lognorm)
  
  h,x,y = hist_and_plot(data[(data['T'] > -20)*(data['T'] < -15)],
                        '-20 < T < -15', yvar, xvar,
                        xlab, ylab, figax=[fig, axs[2,0]],
                        xlim=xlim, ylim=ylim, inverty=False,
                        savename=None, lognorm=lognorm)
  
  h,x,y = hist_and_plot(data[(data['T'] > -25)*(data['T'] < -20)],
                        '-25 < T < -20', yvar, xvar,
                        xlab, ylab, figax=[fig, axs[2,1]],
                        xlim=xlim, ylim=ylim, inverty=False,
                        savename=None, lognorm=lognorm)
  fig.savefig(savename)

title = 'contour_log'
fig, ax = plt.subplots(1,1)
add_one_contour(None, data, [None, None], 'k')
add_one_contour(None, data, [None, 0], 'k')
add_one_contour(ax, data, [0, None], 'C1')
add_one_contour(ax, data, [-5, 0], 'C2')
add_one_contour(ax, data, [-10, -5], 'C3')
add_one_contour(ax, data, [-15, -10], 'C4')
add_one_contour(ax, data, [-20, -15], 'C5')
add_one_contour(ax, data, [-25, -20], 'C6')
add_one_contour(ax, data, [-30, -25], 'C7')

ax.set_xlim([-5,15])
ax.set_ylim([-5,15])
ax.set_title(title)
ax.set_xlabel('DWR Ka W   [dB]')
ax.set_ylabel('DWR X Ka   [dB]')
ax.grid()
ax.legend()
fig.savefig('triple_frequency/' + '_'.join(title.split()) + addlabel + '.png')