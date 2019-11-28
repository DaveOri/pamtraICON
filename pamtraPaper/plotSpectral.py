#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 13:11:37 2019

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

pamX = nc.Dataset('all_hydro/20151119all_hydro_spe_KiXPol.nc')
pamK = nc.Dataset('all_hydro/20151119all_hydro_spe_Joyrad35.nc')
pamW = nc.Dataset('all_hydro/20151119all_hydro_spe_Joyrad94.nc')
pamP = nc.Dataset('20151119hatpro.nc')

rad3 = nc.Dataset('/data/optimice/tripex/tripex_level_02_test/tripex_joy_tricr00_l2_any_v00_20151119000000.nc')
hatp = nc.Dataset('/data/hatpro/jue/data/level1/1511/sups_joy_mwr00_l1_tb_v01_20151119000003.nc')

xfmt = md.DateFormatter('%H')
xDataLim = -1
versus =  1 # Bottom Up

hlim = [0, 8]
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

fig1, axes = plt.subplots(nrows=3, ncols=2, figsize=(10, 13),
                          constrained_layout=True,
                          gridspec_kw={'height_ratios': [3, 1, 1]})
#((ax11, ax12), (ax13, ax14)) = axes
((ax11, ax12), (ax13, ax14), (ax15, ax16)) = axes

tidx = 203#50
hidx = 70#60#18#18#29#70
selTime = tta[tidx, hidx]
selTS = pd.to_datetime(selTime)
selHeight = Ha[tidx, hidx]

print(selTime, selHeight)

rad94file = '/data/hatpro/jue/data/joyrad94/l0/' + selTS.strftime('%Y%m/%d') \
            + '/joyrad94_joyce_' + selTS.strftime('%Y%m%d%H.nc')
radar94 = nc.Dataset(rad94file)
rad94var = radar94.variables
t94 = nc.num2date(rad94var['time'][:], 'seconds since 2001-01-01')
t94idx = np.argmin(np.abs(t94 - selTS))
h94idx = np.argmin(np.abs(rad94var['range'][:] - selHeight))
spec94 = rad94var['spec'][t94idx,:,:]
spec94[spec94 == -999] = np.nan
chirp_idx = list(rad94var['range_offsets'][:] - 1)
chirp_idx.append(spec94.shape[0])
v94 = np.zeros(spec94.shape)
for i,j in enumerate(chirp_idx[:-1]):
  v94[j:chirp_idx[i+1],:] = np.tile(rad94var['velocity'][i,:], [chirp_idx[i+1]-j,1])
r94 = np.tile(rad94var['range'][:][:, np.newaxis], spec94.shape[1])
v94 = -v94

rad35file = '/net/ora/20151119_1633.zspc.nc'
radar35 = nc.Dataset(rad35file)
rad35var = radar35.variables
t35 = nc.num2date(rad35var['time'][:], 'seconds since 1970-01-01 00:00 UTC')
t35idx = np.argmin(np.abs(t35 - selTS))
r35 = rad35var['range'][:]
h35idx = np.argmin(np.abs(r35 - selHeight))

SNRCorFaCo = rad35var['SNRCorFaCo'][t35idx,:]
radConst = rad35var['RadarConst'][t35idx]
npw = rad35var['npw1'][t35idx]
cal = radConst*SNRCorFaCo*(r35/5000.)**2/npw
spec35 = 10.*np.log10(rad35var['SPCco'][t35idx,:,:]*cal[:, np.newaxis])
v35 = sorted(rad35var['doppler'][:])

mesh11 = ax11.pcolormesh(pamK.variables['Radar_Velocity'][:].squeeze(),
                         pamK.variables['height'][tidx,0,:]*0.001,
                         pamK.variables['Radar_Spectrum'][tidx,0,:,0,0,:],
                         vmin=Smin, vmax=Smax, linewidth=0, rasterized=True)
ax11.set_xlim(vlim)
ax11.set_ylim(hlim)
#ax11.set_xlabel('Doppler Velocity   [m/s]')
ax11.get_xaxis().set_ticklabels([])
ax11.set_ylabel('Height   [km]')
ax11.set_title('Model',weight='black')
ax11.text(0.1,0.9,'(a)',transform=ax11.transAxes,weight='black',
                  bbox=dict(facecolor='white'))

mesh12 = ax12.pcolormesh(v35, r35*0.001, spec35, vmin=Smin, vmax=Smax,
                         linewidth=0, rasterized=True)
#mesh12 = ax12.pcolormesh(v94, r94*0.001, spec94)
ax12.set_xlim(vlim)
ax12.set_ylim(hlim)
#ax12.set_xlabel('Doppler Velocity   [m/s]')
ax12.get_xaxis().set_ticklabels([])
ax12.get_yaxis().set_ticklabels([])

ax11.axhline(y=selHeight*0.001, c='r')
ax12.axhline(y=selHeight*0.001, c='r')
ax12.set_title('Observations',weight='black')
ax12.text(0.1,0.9,'(b)',transform=ax12.transAxes,weight='black',
                  bbox=dict(facecolor='white'))
#
#ax13.plot(pamX.variables['Radar_Velocity'][:].squeeze(),
#          pamX.variables['Radar_Spectrum'][tidx,0,hidx,0,0,:],
#          label='X band', c='xkcd:pale red')
ax13.plot(pamK.variables['Radar_Velocity'][:].squeeze(),
          pamK.variables['Radar_Spectrum'][tidx,0,hidx,0,0,:],
          label='Ka band', c='xkcd:medium green')
ax13.plot(pamW.variables['Radar_Velocity'][:].squeeze(),
          pamW.variables['Radar_Spectrum'][tidx,0,hidx,0,0,:],
          label='W band', c='xkcd:denim blue')
linspecW = 10.0**(0.1*pamW.variables['Radar_Spectrum'][tidx,0,hidx,0,0,:])

ax13.set_xlim(vlim)
ax13.set_ylim(splim)
ax13.get_xaxis().set_ticklabels([])
ax13.set_ylabel('Spectral Power   [dB]')
#ax13.set_xlabel('Doppler Velocity   [m/s]')
ax13.text(0.1,0.8,'(c)',transform=ax13.transAxes,weight='black')
#ax13.text(0.5,0.8,'Height= 4.8 km ',transform=ax13.transAxes)
ax13.text(0.5,0.8,'Height= 6 km ',transform=ax13.transAxes)

ax14.plot(v35, spec35[h35idx, :],
          label='Ka band', c='xkcd:medium green')
linspec35 = 10.0**(0.1*spec35[h35idx, :])

ax14.plot(v94[h94idx,:], spec94[h94idx,:],
          label = 'W band', c='xkcd:denim blue')
linspec94 = 10.0**(0.1*spec94[h94idx, :])

ax14.set_xlim(vlim)
ax14.set_ylim(splim)
ax14.get_yaxis().set_ticklabels([])
ax14.get_xaxis().set_ticklabels([])
#ax14.set_xlabel('Doppler Velocity   [m/s]')
ax14.legend()
ax14.text(0.1,0.8,'(d)',transform=ax14.transAxes,weight='black')

hidx = 16#32#18#29#70
selTime = tta[tidx, hidx]
selTS = pd.to_datetime(selTime)
selHeight = Ha[tidx, hidx]
h94idx = np.argmin(np.abs(rad94var['range'][:] - selHeight))
h35idx = np.argmin(np.abs(r35 - selHeight))
print(selTime, selHeight)

#ax15.plot(pamX.variables['Radar_Velocity'][:].squeeze(),
#          pamX.variables['Radar_Spectrum'][tidx,0,hidx,0,0,:],
#          label='X band', c='xkcd:pale red')
ax15.plot(pamK.variables['Radar_Velocity'][:].squeeze(),
          pamK.variables['Radar_Spectrum'][tidx,0,hidx,0,0,:],
          label='Ka band', c='xkcd:medium green')
ax15.plot(pamW.variables['Radar_Velocity'][:].squeeze(),
          pamW.variables['Radar_Spectrum'][tidx,0,hidx,0,0,:],
          label='W band', c='xkcd:denim blue')
linspecW = 10.0**(0.1*pamW.variables['Radar_Spectrum'][tidx,0,hidx,0,0,:])

ax15.set_xlim(vlim)
ax15.set_ylim(splim)
ax15.set_ylabel('Spectral Power   [dB]')
ax15.set_xlabel('Doppler Velocity   [m/s]')
ax15.text(0.1,0.8,'(e)',transform=ax15.transAxes,weight='black')
ax15.text(0.5,0.1,'Height= 0.9 km ',transform=ax15.transAxes)

ax16.plot(v35, spec35[h35idx, :],
          label='Ka band', c='xkcd:medium green')
linspec35 = 10.0**(0.1*spec35[h35idx, :])

ax16.plot(v94[h94idx,:], spec94[h94idx,:],
          label = 'W band', c='xkcd:denim blue')
linspec94 = 10.0**(0.1*spec94[h94idx, :])

ax16.set_xlim(vlim)
ax16.set_ylim(splim)
ax16.get_yaxis().set_ticklabels([])
ax16.set_xlabel('Doppler Velocity   [m/s]')
ax16.legend()
ax16.text(0.1,0.8,'(f)',transform=ax16.transAxes,weight='black')

ax13.grid()
ax14.grid()
ax15.grid()
ax16.grid()

ax11.axhline(y=selHeight*0.001, c='r')
ax12.axhline(y=selHeight*0.001, c='r')
cbar = fig1.colorbar(mesh11, ax=ax12, label='Spectral Power   [dB]')
#fig1.savefig('tripex_spectra_ssrg-rt3.png', dpi=600)
fig1.savefig('tripex_spectra.pdf')
fig1.savefig('tripex_spectra.png', dpi=600)