#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:50:24 2019

@author: dori
"""
import matplotlib.pyplot as plt
plt.close('all')
import matplotlib.dates as md
import pandas as pd
import numpy as np
import netCDF4 as nc

import matplotlib
matplotlib.rcParams.update({'font.size':6})

pamX = nc.Dataset('all_hydro/20151119all_hydro_spe_KiXPol.nc')
pamK = nc.Dataset('all_hydro/20151119all_hydro_spe_Joyrad35.nc')
pamW = nc.Dataset('all_hydro/20151119all_hydro_spe_Joyrad94.nc')
pamP = nc.Dataset('20151119hatpro.nc')

rad3 = nc.Dataset('/data/optimice/tripex/tripex_level_02_test/tripex_joy_tricr00_l2_any_v00_20151119000000.nc')
hatp = nc.Dataset('/data/hatpro/jue/data/level1/1511/sups_joy_mwr00_l1_tb_v01_20151119000003.nc')

xfmt = md.DateFormatter('%m-%d %H')
xfmt = md.DateFormatter('%H')
ylim=(0,12000)
xDataLim = -1
figsize31=(18,18)
figsize21=(18,12)
versus = -1 # Top Down
versus =  1 # Bottom Up

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
  
def hildebrand_sekhon(spectrum, nAvg):
    N = len(spectrum) - 1 # index of the last element
    sortSpec = np.sort(spectrum)
    sortSpec[~np.isfinite(sortSpec)] = 0.0
    sum2 = sortSpec[0:-1].cumsum()**2
    sumSquares= (sortSpec[0:-1]**2).cumsum()
    nPts=np.arange(N)+1.0
    sigma = nAvg*(nPts*sumSquares-sum2)
    a=int(max(nPts[sigma < sum2]))
    
    if(a < 10):
        peak_noise=max(spectrum)
    else:
        peak_noise=sortSpec[a-1]
    signal_detected=1.0*spectrum
    signal_detected[spectrum < peak_noise] = np.nan
    return peak_noise, signal_detected
  
Hx, ttx, Ax, Zex, MDVx, SWx = readPamtra_nc(pamX)
Ha, tta, Aa, Zea, MDVa, SWa = readPamtra_nc(pamK)
Hw, ttw, Aw, Zew, MDVw, SWw = readPamtra_nc(pamW)

fig0 = plt.figure(figsize=(8, 11))
ax1 = plt.subplot2grid((5, 4), (0, 0), colspan=2)
ax2 = plt.subplot2grid((5, 4), (1, 0), colspan=2)
ax3 = plt.subplot2grid((5, 4), (2, 0), colspan=2)
ax4 = plt.subplot2grid((5, 4), (3, 0), colspan=2)
ax5 = plt.subplot2grid((5, 4), (4, 0), colspan=2)
ax6 = plt.subplot2grid((5, 4), (0, 2), colspan=2)
ax7 = plt.subplot2grid((5, 4), (1, 2), colspan=2)
ax8 = plt.subplot2grid((5, 4), (2, 2), colspan=2)
ax9 = plt.subplot2grid((5, 4), (3, 2), colspan=2)
ax10 = plt.subplot2grid((5, 4), (4, 2), colspan=2)

fig1 = plt.figure(figsize=(8, 8))
ax11 = plt.subplot2grid((5, 2), (0, 0), rowspan=4)
ax12 = plt.subplot2grid((5, 2), (0, 1), rowspan=4)
ax13 = plt.subplot2grid((5, 2), (4, 0), colspan=1)
ax14 = plt.subplot2grid((5, 2), (4, 1), colspan=1)
#axCB = plt.subplot2grid((5, 3), (0, 2), colspan=3)

hlim = [0, 10]
TlimK = [0, 150]
TlimV = [100, 300]
vlim = [-2, 8]
splim = [-70, 0]
Smin, Smax = splim
Zmin, Zmax = -35, 25
Vmin, Vmax = -1, 5
Dmin, Dmax = -5, 20

mesh1 = ax1.pcolormesh(getTime(pamK, 'datatime'),
                       pamK.variables['height'][0,0,:]*0.001,
                       (Zea-Aa).T, cmap='jet', vmin=Zmin, vmax=Zmax)
ax1.set_ylim(hlim)
ax1.get_xaxis().set_ticks([])
ax1.set_ylabel('Height    [km]')

mesh2 = ax2.pcolormesh(getTime(pamK, 'datatime'),
                       pamK.variables['height'][0,0,:]*0.001,
                       -MDVa.T, cmap='RdBu', vmin=Vmin, vmax=Vmax)
ax2.set_ylim(hlim)
ax2.get_xaxis().set_ticks([])
ax2.set_ylabel('Height    [km]')

mesh3 = ax3.pcolormesh(getTime(pamK, 'datatime'),
                       pamK.variables['height'][0,0,:]*0.001,
                       (Zea-Aa-Zew+Aw).T, cmap='nipy_spectral', vmin=Dmin, vmax=Dmax)
ax3.set_ylim(hlim)
ax3.set_ylabel('Height    [km]')
ax3.set_xlabel('time')
ax3.xaxis.set_major_formatter(xfmt)


ax4.plot(getTime(pamP, 'datatime'), pamP.variables['tb'][:,0,1,31,:7,0]) # downwelling at 0 meters)
ax4.legend(f2labels(pamP.variables['frequency'][:7]))
ax4.set_ylim(TlimK)
ax4.set_xlim(ax1.get_xlim())
ax4.get_xaxis().set_ticks([])
ax4.set_ylabel('T$_b$   [K]')

ax5.plot(getTime(pamP, 'datatime'), pamP.variables['tb'][:,0,1,31,7:,0]) # downwelling at 0 meters)
ax5.legend(f2labels(pamP.variables['frequency'][7:]))
ax5.set_ylim(TlimV)
ax5.set_xlim(ax1.get_xlim())
ax5.set_ylabel('T$_b$   [K]')
ax5.set_xlabel('time')
ax5.xaxis.set_major_formatter(xfmt)

mesh6 = ax6.pcolormesh(getTime(rad3, 'time'),
                       rad3.variables['range'][:]*0.001,
                       rad3.variables['dbz_ka'][:].T, cmap='jet',
                       vmin=Zmin, vmax=Zmax)
ax6.set_ylim(hlim)
ax6.get_xaxis().set_ticks([])
ax6.get_yaxis().set_ticks([])


mesh7 = ax7.pcolormesh(getTime(rad3, 'time'),
                       rad3.variables['range'][:]*0.001,
                       -rad3.variables['rv_ka'][:].T, cmap='RdBu',
                       vmin=Vmin, vmax=Vmax)
ax7.set_ylim(hlim)
ax7.get_xaxis().set_ticks([])
ax7.get_yaxis().set_ticks([])

mesh8 = ax8.pcolormesh(getTime(rad3, 'time'),
                       rad3.variables['range'][:]*0.001,
                       rad3.variables['dbz_ka'][:].T-rad3.variables['dbz_w'][:].T,
                       cmap='nipy_spectral', vmin=Dmin, vmax=Dmax)
ax8.set_ylim(hlim)
ax8.set_xlim(ax1.get_xlim())
ax8.get_yaxis().set_ticks([])
ax8.set_xlabel('time')
ax8.xaxis.set_major_formatter(xfmt)

ele = hatp.variables['ele'][:]
elemask = np.abs(ele-90.) < 1.
ax9.plot(getTime(hatp, 'time')[elemask], hatp.variables['tb'][elemask, :7])
ax9.set_ylim(TlimK)
ax9.set_xlim(ax1.get_xlim())
ax9.get_xaxis().set_ticks([])
ax9.get_yaxis().set_ticks([])

ax10.plot(getTime(hatp, 'time')[elemask], hatp.variables['tb'][elemask, 7:])
ax10.set_ylim(TlimV)
ax10.set_xlim(ax1.get_xlim())
ax10.get_yaxis().set_ticks([])
ax10.set_xlabel('time')
ax10.xaxis.set_major_formatter(xfmt)

fig0.subplots_adjust(right=0.7)
cax6 = fig0.add_axes([0.85, 0.75, 0.05, 0.15])
cax7 = fig0.add_axes([0.85, 0.55, 0.05, 0.15])
cax8 = fig0.add_axes([0.85, 0.35, 0.05, 0.15])
plt.colorbar(mesh6, cax=cax6, label='Z$_e$   [dBZ]')
plt.colorbar(mesh7, cax=cax7, label='Mean Doppler Velocity  [m/s]')
plt.colorbar(mesh8, cax=cax8, label='DWR$_{K_aW}$  [dB]')

#fig0.tight_layout()
fig0.savefig('tripex_plots.png', dpi=600)

tidx = 197#939#6572#4800
hidx = 42#18#18#29#70
selTime = tta[tidx, hidx]
selTS = pd.to_datetime(selTime)
selHeight = Ha[tidx, hidx]

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
                         vmin=Smin, vmax=Smax)
ax11.set_xlim(vlim)
ax11.set_ylim(hlim)
#ax11.set_xlabel('Doppler Velocity   [m/s]')
ax11.get_xaxis().set_ticks([])
ax11.set_ylabel('Height   [km]')

mesh12 = ax12.pcolormesh(v35, r35*0.001, spec35, vmin=Smin, vmax=Smax)
#mesh12 = ax12.pcolormesh(v94, r94*0.001, spec94)
ax12.set_xlim(vlim)
ax12.set_ylim(hlim)
#ax12.set_xlabel('Doppler Velocity   [m/s]')
ax12.get_xaxis().set_ticks([])
ax12.get_yaxis().set_ticks([])
#plt.colorbar(mesh12, ax=ax12, label='Spectral Power   [dB]')

ax13.plot(pamX.variables['Radar_Velocity'][:].squeeze(),
          pamX.variables['Radar_Spectrum'][tidx,0,hidx,0,0,:],
          label='X band', c='xkcd:pale red')
ax13.plot(pamK.variables['Radar_Velocity'][:].squeeze(),
          pamK.variables['Radar_Spectrum'][tidx,0,hidx,0,0,:],
          label='Ka band', c='xkcd:medium green')
ax13.plot(pamW.variables['Radar_Velocity'][:].squeeze(),
          pamW.variables['Radar_Spectrum'][tidx,0,hidx,0,0,:],
          label='W band', c='xkcd:denim blue')
linspecW = 10.0**(0.1*pamW.variables['Radar_Spectrum'][tidx,0,hidx,0,0,:])
peaknoiseW, detectedSpectrumW = hildebrand_sekhon(linspecW, 20)
ax13.legend()
ax13.set_xlim(vlim)
ax13.set_ylim(splim)
ax13.set_ylabel('Spectral Power   [dB]')
ax13.set_xlabel('Doppler Velocity   [m/s]')

ax14.plot(v35, spec35[h35idx, :],
          label='Ka band', c='xkcd:medium green')
linspec35 = 10.0**(0.1*spec35[h35idx, :])
peaknoise35, detectedSpectrum35 = hildebrand_sekhon(linspec35, 20)

ax14.plot(v94[h94idx,:], spec94[h94idx,:],
          label = 'W band', c='xkcd:denim blue')
linspec94 = 10.0**(0.1*spec94[h94idx, :])
peaknoise94, detectedSpectrum94 = hildebrand_sekhon(linspec94, 17)

ax14.set_xlim(vlim)
ax14.set_ylim(splim)
#ax14.set_ylabel('Spectral Power   [dB]')
ax14.get_yaxis().set_ticks([])
ax14.set_xlabel('Doppler Velocity   [m/s]')

fig1.subplots_adjust(bottom=0.5, left=0.15, right=0.7)
axCB = fig1.add_axes([0.85, 0.3, 0.05, 0.5])
cbar = fig1.colorbar(mesh11, cax=axCB, label='Spectral Power   [dB]')
#fig0.tight_layout()

print('X:', Zex[tidx, hidx], '  ',10.0*np.log10((10.0**(0.1*pamX.variables['Radar_Spectrum'][tidx,0,hidx,0,0,:])).sum()))
print('K:', Zea[tidx, hidx], '  ',10.0*np.log10((10.0**(0.1*pamK.variables['Radar_Spectrum'][tidx,0,hidx,0,0,:])).sum()))
print('W:', Zew[tidx, hidx], '  ',10.0*np.log10((10.0**(0.1*pamW.variables['Radar_Spectrum'][tidx,0,hidx,0,0,:])).sum()))

print('K:',10.0*np.log10(rad35var['Zg'][t35idx, h35idx]), '     ',10.0*np.log10((10**(0.1*spec35[h35idx,:])).sum()))
print('W:',10.0*np.log10(rad94var['Ze'][t94idx, h94idx]), '     ',10.0*np.log10((10**(0.1*spec94[h94idx,:])).sum()))

print('Noise:', 10.0*np.log10(peaknoise35), 10.0*np.log10(peaknoise94) )
fig1.savefig('tripex_spectra.png', dpi=600)

#noise_profile = 0.0*SNRCorFaCo
#for i in range(spec35.shape[0]):
#  linspec35 = 10.0**(0.1*spec35[i, :])
#  peaknoise35, detectedSpectrum35 = hildebrand_sekhon(linspec35, 20)
#  noise_profile[i] = 10.0*np.log10(peaknoise35)
#
#def noiseFunc(r, N0):
#  return N0 + 20*np.log10(r/1000.)
#
#from scipy.optimize import curve_fit
#mask = (0.0*r35.data).astype(bool)
#mask[:] = True
#mask[:10] = False
#mask[143:231] = False
#pars, covs = curve_fit(noiseFunc, r35[mask], noise_profile[mask])
#plt.figure()
#plt.scatter(r35[mask]/1000., noise_profile[mask])
#plt.scatter(r35[~mask]/1000., noise_profile[~mask])
#plt.plot(r35/1000., noiseFunc(r35, *pars))
#
#noise_profile = 0.0*r94[:,0]
#for i in range(spec94.shape[0]):
#  linspec94 = 10.0**(0.1*spec94[i, :])
#  try:
#    peaknoise94, detectedSpectrum94 = hildebrand_sekhon(linspec94, 17)
#    noise_profile[i] = 10.0*np.log10(peaknoise94)
#  except:
#    noise_profile[i] = np.nan
#
#mask = (0.0*r94.data[:,0]).astype(bool)
#mask[:] = np.isfinite(noise_profile)
##mask[:10] = False
##mask[143:231] = False
#pars, covs = curve_fit(noiseFunc, r94[mask,0], noise_profile[mask])
#plt.figure()
#plt.scatter(r94[mask,0]/1000., noise_profile[mask])
#plt.scatter(r94[~mask,0]/1000., noise_profile[~mask])
#plt.plot(r94[:,0]/1000., noiseFunc(r94[:,0], *pars))