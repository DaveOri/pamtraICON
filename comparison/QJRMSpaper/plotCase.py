#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 12:13:52 2019

@author: dori
"""

import matplotlib.pyplot as plt
plt.close('all')
import matplotlib.dates as md
import pandas as pd
import numpy as np
import netCDF4 as nc

date = '20151119'
oIce = '/data/optimice/'
radPath = oIce + 'tripex/tripex_level_02_test/tripex_joy_tricr00_l2_any_v00_'
rad3 = nc.Dataset(radPath + date + '000000.nc')

pamPath = oIce + '/pamtra_runs/tripex/data/all_hydro/'
pamX = nc.Dataset(pamPath + date + 'all_hydro_mom_KiXPol.nc')
pamK = nc.Dataset(pamPath + date + 'all_hydro_mom_Joyrad35.nc')
pamW = nc.Dataset(pamPath + date + 'all_hydro_mom_Joyrad94.nc')

plt.rcParams.update({'font.size': 9})

xfmt = md.DateFormatter('%H')
xDataLim = -1
versus =  1 # Bottom Up

hlim = [0, 10]
TlimK = [0, 120]
TlimV = [100, 300]
vlim = [-2, 8]
splim = [-70, 0]
Smin, Smax = splim
Zmin, Zmax = -35, 25
Vmin, Vmax = -1, 5
Dmin, Dmax = -5, 20

# Extract netCDF4 time variables
def getTime(ncdata, timevar):
  return nc.num2date(times=ncdata.variables[timevar][:],
                     units=ncdata.variables[timevar].units).squeeze()

def f2labels(frequencies):
  return [str(f)+' GHz' for f in frequencies]

# Define Plotting Function
def plot_variable(x,y,v,axes,
                  xlab=None,ylab=None,vlab=None,title=None,
                  vmin=None,vmax=None,xlim=None,ylim=None,
                  cmap='jet'):
    mesh = axes.pcolormesh(x, y, v, vmin=vmin, vmax=vmax, cmap=cmap,
                           linewidth=0, rasterized=True)
    if title is not None:
        axes.text(0.1,0.9,title,transform=axes.transAxes,weight='black',
                  bbox=dict(facecolor='white'))
    plt.colorbar(mesh,label=vlab,ax=axes)
    if xlab is not None:
        axes.set_xlabel(xlab)
    if ylab is not None:
        axes.set_ylabel(ylab)
    axes.set_xlim(xlim)
    axes.set_ylim(ylim)

def readPamtra_nc(ncfile):
    #runDataset = Dataset(ncfile)
    runVars = ncfile.variables
    H = (runVars['height'][:,0,:])[:xDataLim,:]
    ttt = pd.to_datetime(runVars['datatime'][:,0],unit='s')
    tt = (np.tile(ttt,(H.shape[1],1)).T)[:xDataLim,:]
    print(tt.shape, H.shape)
    a = 2.0*(runVars['Attenuation_Hydrometeors'][:,0,:,0,0] + runVars['Attenuation_Atmosphere'][:,0,:,0])
    A = a[:,::versus].cumsum(axis=1)[:,::versus][:xDataLim,:]
    Ze = runVars['Ze'][:,0,:,0,0,0][:xDataLim,:]
    MDV = -runVars['Radar_MeanDopplerVel'][:,0,:,0,0,0][:xDataLim,:]
    SW = runVars['Radar_SpectrumWidth'][:,0,:,0,0,0][:xDataLim,:]
    return H, tt, A, Ze, MDV, SW

Hx, ttx, Ax, Zex, MDVx, SWx = readPamtra_nc(pamX)
Ha, tta, Aa, Zea, MDVa, SWa = readPamtra_nc(pamK)
Hw, ttw, Aw, Zew, MDVw, SWw = readPamtra_nc(pamW)

fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(9.5, 6),
                         sharex=True,
                         constrained_layout=True)
((ax0, ax1), (ax2, ax3), (ax4, ax5)) = axes

mesh0 = ax0.pcolormesh(getTime(pamX, 'datatime'),
                       pamX.variables['height'][0,0,:]*0.001,
                       (Zex-Ax).T, cmap='jet', vmin=Zmin, vmax=Zmax,
                       linewidth=0, rasterized=True)
ax0.set_ylim(hlim)
#ax1.get_xaxis().set_ticklabels([])
ax0.set_ylabel('Height    [km]')
#ax0.xaxis.set_major_formatter(xfmt)
ax0.grid()
ax0.text(0.1, 0.8, 'Simulated X band',
         transform=ax0.transAxes, weight='normal',
         bbox=dict(facecolor='white'))

mesh1 = ax1.pcolormesh(getTime(rad3, 'time'),
                       rad3.variables['range'][:]*0.001,
                       rad3.variables['dbz_x'][:].T, cmap='jet',
                       vmin=Zmin, vmax=Zmax,
                       linewidth=0, rasterized=True)
ax1.set_ylim(hlim)
#ax1.get_xaxis().set_ticklabels([])
ax1.get_yaxis().set_ticklabels([])
#ax1.xaxis.set_major_formatter(xfmt)
ax1.grid()
ax1.text(0.1, 0.8, 'KiXPol',
         transform=ax1.transAxes, weight='normal',
         bbox=dict(facecolor='white'))

mesh2 = ax2.pcolormesh(getTime(pamK, 'datatime'),
                       pamK.variables['height'][0,0,:]*0.001,
                       (Zea-Aa).T, cmap='jet', vmin=Zmin, vmax=Zmax,
                       linewidth=0, rasterized=True)
ax2.set_ylim(hlim)
#ax2.get_xaxis().set_ticklabels([])
ax2.set_ylabel('Height    [km]')
#ax2.xaxis.set_major_formatter(xfmt)
ax2.grid()
ax2.text(0.1, 0.8, 'Simulated Ka band',
         transform=ax2.transAxes, weight='normal',
         bbox=dict(facecolor='white'))

mesh3 = ax3.pcolormesh(getTime(rad3, 'time'),
                       rad3.variables['range'][:]*0.001,
                       rad3.variables['dbz_ka'][:].T, cmap='jet',
                       vmin=Zmin, vmax=Zmax,
                       linewidth=0, rasterized=True)
ax3.set_ylim(hlim)
#ax3.get_xaxis().set_ticklabels([])
ax3.get_yaxis().set_ticklabels([])
#ax3.xaxis.set_major_formatter(xfmt)
ax3.grid()
ax3.text(0.1, 0.8, 'Joyrad35',
         transform=ax3.transAxes, weight='normal',
         bbox=dict(facecolor='white'))


mesh4 = ax4.pcolormesh(getTime(pamW, 'datatime'),
                       pamW.variables['height'][0,0,:]*0.001,
                       (Zew-Aw).T, cmap='jet', vmin=Zmin, vmax=Zmax,
                       linewidth=0, rasterized=True)
ax4.set_ylim(hlim)
ax4.set_ylabel('Height    [km]')
ax4.set_xlabel('time')
ax4.xaxis.set_major_formatter(xfmt)
ax4.grid()
ax4.text(0.1, 0.8, 'Simulated W band',
         transform=ax4.transAxes, weight='normal',
         bbox=dict(facecolor='white'))


mesh5 = ax5.pcolormesh(getTime(rad3, 'time'),
                       rad3.variables['range'][:]*0.001,
                       rad3.variables['dbz_w'][:].T,
                       cmap='jet', vmin=Zmin, vmax=Zmax,
                       linewidth=0, rasterized=True)
ax5.set_ylim(hlim)
ax5.set_xlim(ax0.get_xlim())
ax5.get_yaxis().set_ticklabels([])
ax5.set_xlabel('time')
ax5.xaxis.set_major_formatter(xfmt)
ax5.grid()
ax5.text(0.1, 0.8, 'Joyrad94',
         transform=ax5.transAxes, weight='normal',
         bbox=dict(facecolor='white'))


plt.colorbar(mesh1, ax=ax1, label='Z$_e$   [dBZ]',
             ticks=[-35, -20, -5, 10, 25])
plt.colorbar(mesh3, ax=ax3, label='Z$_e$   [dBZ]',
             ticks=[-35, -20, -5, 10, 25])
plt.colorbar(mesh5, ax=ax5, label='Z$_e$   [dBZ]',
             ticks=[-35, -20, -5, 10, 25])
fig.suptitle('Radar reflectivity 2015-11-19',
             fontsize=12, fontweight='heavy')
#fig.tight_layout(pad=1.5, h_pad=0.5, w_pad=0.5)
fig.savefig(date + 'tripex_plots.pdf')
fig.savefig(date + 'tripex_plots.png', dpi=350)

