#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 13:10:07 2019

@author: dori
"""

import matplotlib.pyplot as plt
plt.close('all')
import matplotlib.dates as md
import pandas as pd
import numpy as np
import netCDF4 as nc

import matplotlib
matplotlib.rcParams.update({'font.size':14})
matplotlib.rcParams.update({'legend.fontsize':12})

pamX = nc.Dataset('all_hydro/20151119all_hydro_spe_KiXPol.nc')
pamK = nc.Dataset('all_hydro/20151119all_hydro_spe_Joyrad35.nc')
pamW = nc.Dataset('all_hydro/20151119all_hydro_spe_Joyrad94.nc')
pamP = nc.Dataset('20151119hatpro.nc')

rad3 = nc.Dataset('/data/optimice/tripex/tripex_level_02_test/tripex_joy_tricr00_l2_any_v00_20151119000000.nc')
hatp = nc.Dataset('/data/hatpro/jue/data/level1/1511/sups_joy_mwr00_l1_tb_v01_20151119000003.nc')

flag = hatp.variables['flag'][:]
flags = ['{0:010b}'.format(i) for i in flag.data]
rainflag = [bool(int(i[-4])) for i in flags]

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
    mesh = axes.pcolormesh(x,y,v,vmin=vmin,vmax=vmax,cmap=cmap)
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

fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(10, 13), sharex=True,
                         constrained_layout=True)
((ax0, ax1), (ax2, ax3), (ax4, ax5), (ax6, ax7), (ax8, ax9)) = axes

mesh0 = ax0.pcolormesh(getTime(pamK, 'datatime'),
                       pamK.variables['height'][0,0,:]*0.001,
                       (Zea-Aa).T, cmap='jet', vmin=Zmin, vmax=Zmax)
ax0.set_ylim(hlim)
#ax1.get_xaxis().set_ticklabels([])
ax0.set_ylabel('Height    [km]')

mesh1 = ax1.pcolormesh(getTime(rad3, 'time'),
                       rad3.variables['range'][:]*0.001,
                       rad3.variables['dbz_ka'][:].T, cmap='jet',
                       vmin=Zmin, vmax=Zmax)
ax1.set_ylim(hlim)
#ax1.get_xaxis().set_ticklabels([])
ax1.get_yaxis().set_ticklabels([])

mesh2 = ax2.pcolormesh(getTime(pamK, 'datatime'),
                       pamK.variables['height'][0,0,:]*0.001,
                       -MDVa.T, cmap='RdBu', vmin=Vmin, vmax=Vmax)
ax2.set_ylim(hlim)
#ax2.get_xaxis().set_ticklabels([])
ax2.set_ylabel('Height    [km]')

mesh3 = ax3.pcolormesh(getTime(rad3, 'time'),
                       rad3.variables['range'][:]*0.001,
                       -rad3.variables['rv_ka'][:].T, cmap='RdBu',
                       vmin=Vmin, vmax=Vmax)
ax3.set_ylim(hlim)
#ax3.get_xaxis().set_ticklabels([])
ax3.get_yaxis().set_ticklabels([])

mesh4 = ax4.pcolormesh(getTime(pamK, 'datatime'),
                       pamK.variables['height'][0,0,:]*0.001,
                       (Zea-Aa-Zew+Aw).T, cmap='nipy_spectral', vmin=Dmin, vmax=Dmax)
ax4.set_ylim(hlim)
ax4.set_ylabel('Height    [km]')
ax4.set_xlabel('time')
ax4.xaxis.set_major_formatter(xfmt)

mesh5 = ax5.pcolormesh(getTime(rad3, 'time'),
                       rad3.variables['range'][:]*0.001,
                       rad3.variables['dbz_ka'][:].T-rad3.variables['dbz_w'][:].T,
                       cmap='nipy_spectral', vmin=Dmin, vmax=Dmax)
ax5.set_ylim(hlim)
ax5.set_xlim(ax0.get_xlim())
ax5.get_yaxis().set_ticklabels([])
ax5.set_xlabel('time')
ax5.xaxis.set_major_formatter(xfmt)

## PASSIVE

ax6.plot(getTime(pamP, 'datatime'), pamP.variables['tb'][:,0,1,31,:7,0]) # downwelling at 0 meters)
ax6.legend(f2labels(pamP.variables['frequency'][:7]))
ax6.set_ylim(TlimK)
ax6.set_xlim(ax0.get_xlim())
#ax6.get_xaxis().set_ticklabels([])
ax6.set_ylabel('T$_b$   [K]')

ele = hatp.variables['ele'][:]
elemask = np.abs(ele-90.) < 1.
hatprotime = getTime(hatp, 'time')
rain = 5.0*np.ones(hatprotime.shape)
rainmasked = np.ma.MaskedArray(data = rain,
                               mask = ~np.array(rainflag),
                               fill_value = np.nan)
ax7.plot(hatprotime[elemask], hatp.variables['tb'][elemask, :7])
ax7.plot(hatprotime, rainmasked, c='k', lw=3, label='rain flag')
ax7.set_ylim(TlimK)
ax7.set_xlim(ax0.get_xlim())
#ax7.get_xaxis().set_ticklabels([])
ax7.get_yaxis().set_ticklabels([])
ax7.legend(loc=2)

ax8.plot(getTime(pamP, 'datatime'), pamP.variables['tb'][:,0,1,31,7:,0]) # downwelling at 0 meters)
ax8.legend(f2labels(pamP.variables['frequency'][7:]))
ax8.set_ylim(TlimV)
ax8.set_xlim(ax0.get_xlim())
ax8.set_ylabel('T$_b$   [K]')
ax8.set_xlabel('time')
ax8.xaxis.set_major_formatter(xfmt)

ax9.plot(getTime(hatp, 'time')[elemask], hatp.variables['tb'][elemask, 7:])
ax9.plot(hatprotime, 100+rainmasked, c='k', lw=3, label='rain flag')
ax9.set_ylim(TlimV)
ax9.set_xlim(ax1.get_xlim())
ax9.get_yaxis().set_ticks([])
ax9.set_xlabel('time')
ax9.xaxis.set_major_formatter(xfmt)

plt.colorbar(mesh1, ax=ax1, label='Z$_e$   [dBZ]')
plt.colorbar(mesh3, ax=ax3, label='MDV  [m/s]')
plt.colorbar(mesh5, ax=ax5, label='DWR$_{K_aW}$  [dB]')
fig.savefig('tripex_plots.png', dpi=600)