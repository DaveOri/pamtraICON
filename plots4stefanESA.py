# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 17:34:01 2017

@author: dori
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.dates as md
import numpy as np
import pandas as pd

plt.close('all')

# Define some paths and filenames

icon_folder = '/data/inscape/icon/experiments/'

mid_lat_filename   = icon_folder + 'tripex_220km/newicon/METEOGRAM_patch001_joyce.nc'
mid_lat_res_filename = 'ICON/newICON5f20151124_dfNEW.nc'
mix_phase_filename = icon_folder + 'acloud/RF06-2705/postproc/rf06_sat_all.nc'
mix_phase_res_filename = 'ICON4vera/mixed_phase_5f.nc'

# Open the netcdf files

mid_lat_data = Dataset(mid_lat_filename)
mid_lat_res = Dataset(mid_lat_res_filename)
mix_phase_data = Dataset(mix_phase_filename)
mix_phase_res = Dataset(mix_phase_res_filename)

# Extract Variables and Dimensions

ml_dt_dim = mid_lat_data.dimensions
ml_dt_var = mid_lat_data.variables
ml_rs_dim = mid_lat_res.dimensions
ml_rs_var = mid_lat_res.variables

mp_dt_dim = mix_phase_data.dimensions
mp_dt_var = mix_phase_data.variables
mp_rs_dim = mix_phase_res.dimensions
mp_rs_var = mix_phase_res.variables

# Define Plotting Function

def plot_variable(x,y,v,axes,
                  xlab=None,ylab=None,vlab=None,title=None,
                  vmin=None,vmax=None,xlim=None,ylim=None):
    mesh = axes.pcolormesh(x,y,v,vmin=vmin,vmax=vmax,cmap='jet')
    axes.set_title(title)
    plt.colorbar(mesh,label=vlab,ax=axes)
    if xlab is not None:
        axes.set_xlabel(xlab)
    if ylab is not None:
        axes.set_ylabel(ylab)
    axes.set_xlim(xlim)
    axes.set_ylim(ylim)

###############################################################################
# Plot Z

H = ml_rs_var['height'][:,0,:]*0.001
times   = ml_dt_var['time'] # seconds since 2015-11-24 02:00:03 proplectic gregorian
units=times.units.split('since')[0]
basetime=pd.to_datetime(times.units.split('since')[-1])
dtimes = pd.to_timedelta(times[:],unit=str(units)[0])
tt=np.tile((basetime + dtimes),(H.shape[1],1)).T

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
ylim=(0,15)
xfmt = md.DateFormatter('%H:%M')
Zx = ml_rs_var['Ze'][:,0,:,0,0,0]
Zu = ml_rs_var['Ze'][:,0,:,1,0,0]
Za = ml_rs_var['Ze'][:,0,:,2,0,0]
Zw = ml_rs_var['Ze'][:,0,:,3,0,0]
Zg = ml_rs_var['Ze'][:,0,:,4,0,0]
plot_variable(tt,H,Zu,ax1,None,'height [km]','dBZ','Ku-band',-20,20,ylim=ylim)
plot_variable(tt,H,Za,ax2,None,'height [km]','dBZ','Ka-band',-20,20,ylim=ylim)
plot_variable(tt,H,Zw,ax3,None,'height [km]','dBZ', 'W-band',-20,20,ylim=ylim)
plot_variable(tt,H,Zg,ax4,'time','height [km]','dBZ', 'G-band',-20,20,ylim=ylim)

ax4.xaxis.set_major_formatter(xfmt)
f.tight_layout(pad=0)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
DWRua = Zu-Za
DWRaw = Za-Zw
DWRag = Za-Zg
DWRwg = Zw-Zg
plot_variable(tt,H,DWRua,ax1,None,'height [km]','dB','DWR$_{Ku Ka}$',0,20,ylim=ylim)
plot_variable(tt,H,DWRaw,ax2,None,'height [km]','dB','DWR$_{Ka W}$',0,20,ylim=ylim)
plot_variable(tt,H,DWRag,ax3,None,'height [km]','dB','DWR$_{Ka G}$',0,20,ylim=ylim)
plot_variable(tt,H,DWRwg,ax4,'time','height [km]','dB','DWR$_{W G}$',0,20,ylim=ylim)
ax4.xaxis.set_major_formatter(xfmt)
f.tight_layout(pad=0)


###############################################################################
# Plot hydrometeors

f,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,sharey=True)
var = ml_dt_var['QNI']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
var_dim_labels = var.dimensions
xvar = ml_dt_var[var_dim_labels[0]]
yvar = ml_dt_var[var_dim_labels[1]]
xval, yval = np.meshgrid(xvar[:], 0.001*yvar[:])
plot_variable(tt,H,vval,ax1,None,'height [Km]',var.name+'  '+var.unit,var.long_name)
var = ml_dt_var['QNS']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
plot_variable(tt,H,vval,ax2,None,None,var.name+'  '+var.unit,var.long_name)
var = ml_dt_var['QI']
vval = np.ma.masked_less(np.flip(var[:],1),1e-7)
plot_variable(tt,H,vval,ax3,'time','height [Km]',var.name+'  '+var.unit,var.long_name)
var = ml_dt_var['QS']
vval = np.ma.masked_less(np.flip(var[:],1),1e-7)
plot_variable(tt,H,vval,ax4,'time',None,var.name+'  '+var.unit,var.long_name)
ax4.xaxis.set_major_formatter(xfmt)

###############################################################################
##      MIXED PHASE CASE
###############################################################################
# Plot hydrometeors

var = mp_dt_var['QNI']
vval = np.ma.masked_less(var[:],0.1) #var[:].T 
var_dim_labels = var.dimensions
xvar = mp_dt_var[var_dim_labels[1]]
yvar = mp_dt_var['height_m']
xval, yval = np.meshgrid(xvar[:], 0.001*yvar[:])

f,((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2)
ylim = (0,5)
plot_variable(xval,yval,vval,ax1,None,'height [Km]','QNI kg-1','number concentration ice',ylim=ylim)
var = mp_dt_var['QNS']
vval = np.ma.masked_less(var[:],0.1) #var[:].T
plot_variable(xval,yval,vval,ax3,None,'height [Km]','QNS kg-1','number concentration snow',ylim=ylim)
var = mp_dt_var['QNC']
vval = np.ma.masked_less(var[:],0.1) #var[:].T
plot_variable(xval,yval,vval,ax5,xvar.name,'height [Km]','QNC kg-1','number concentration cloud droplets',ylim=ylim)

var = mp_dt_var['QI']
vval = np.ma.masked_less(var[:],1e-8) #var[:].T
plot_variable(xval,yval,vval,ax2,None,'height [Km]','QI  kg kg-1','specific ice content' ,ylim=ylim)
var = mp_dt_var['QS']
vval = np.ma.masked_less(var[:],1e-8) #var[:].T
plot_variable(xval,yval,vval,ax4,None,'height [Km]','QS  kg kg-1','snow mixing ratio' ,ylim=ylim)
var = mp_dt_var['QC']
vval = np.ma.masked_less(var[:],1e-8) #var[:].T
plot_variable(xval,yval,vval,ax6,xvar.name,'height [Km]','QC  kg kg-1','cloud droplets mixing ratio',ylim=ylim)

###############################################################################
# Plot Z

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
Zx = mp_rs_var['Ze'][:,0,:,0,0,0]
Zu = mp_rs_var['Ze'][:,0,:,1,0,0]
Za = mp_rs_var['Ze'][:,0,:,2,0,0]
Zw = mp_rs_var['Ze'][:,0,:,3,0,0]
Zg = mp_rs_var['Ze'][:,0,:,4,0,0]
plot_variable(xval,yval,np.flip(Zu.T,0),ax1,None,'height [km]','dBZ','Ku-band',-20,20,ylim=ylim)
plot_variable(xval,yval,np.flip(Za.T,0),ax2,None,'height [km]','dBZ','Ka-band',-20,20,ylim=ylim)
plot_variable(xval,yval,np.flip(Zw.T,0),ax3,None,'height [km]','dBZ', 'W-band',-20,20,ylim=ylim)
plot_variable(xval,yval,np.flip(Zg.T,0),ax4,'lat','height [km]','dBZ', 'G-band',-20,20,ylim=ylim)
f.tight_layout(pad=0)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
DWRua = Zu-Za
DWRaw = Za-Zw
DWRag = Za-Zg
DWRwg = Zw-Zg
plot_variable(xval,yval,np.flip(DWRua.T,0),ax1,None,'height [km]','dBZ','DWR$_{Ku Ka}$',0,20,ylim=ylim)
plot_variable(xval,yval,np.flip(DWRaw.T,0),ax2,None,'height [km]','dBZ','DWR$_{Ka W}$',0,20,ylim=ylim)
plot_variable(xval,yval,np.flip(DWRag.T,0),ax3,None,'height [km]','dBZ', 'DWR$_{Ka G}$',0,20,ylim=ylim)
plot_variable(xval,yval,np.flip(DWRwg.T,0),ax4,'lat','height [km]','dBZ', 'DWR$_{W G}$',0,20,ylim=ylim)
f.tight_layout(pad=0)
