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
nya_filename = icon_folder + 'nyalesund/newicon-2017-06-23-albedo/METEOGRAM_patch004_awipev.nc'
nya_res_filename = 'ICON4vera/nyalesund_patch004_5f.nc'

frt_filename = icon_folder + 'fronts_postproc/METEOGRAM_patch004_joyce_26only.nc'
frt_res_filename = 'ICON4vera/fronts26only.nc'


# Open the netcdf files

mid_lat_data = Dataset(mid_lat_filename)
mid_lat_res = Dataset(mid_lat_res_filename)
mix_phase_data = Dataset(mix_phase_filename)
mix_phase_res = Dataset(mix_phase_res_filename)
nya_data = Dataset(nya_filename)
nya_res = Dataset(nya_res_filename)
frt_data = Dataset(frt_filename)
frt_res = Dataset(frt_res_filename)

# Extract Variables and Dimensions

ml_dt_dim = mid_lat_data.dimensions
ml_dt_var = mid_lat_data.variables
ml_rs_dim = mid_lat_res.dimensions
ml_rs_var = mid_lat_res.variables

ny_dt_dim = nya_data.dimensions
ny_dt_var = nya_data.variables
ny_rs_dim = nya_res.dimensions
ny_rs_var = nya_res.variables

fr_dt_dim = frt_data.dimensions
fr_dt_var = frt_data.variables
fr_rs_dim = frt_res.dimensions
fr_rs_var = frt_res.variables

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

versus = -1 # Top Down
#versus =  1 # Bottom Up

###############################################################################
# Plot Z

H = fr_rs_var['height'][:,0,:]*0.001
times   = fr_dt_var['time'] # seconds since 2015-06-23 00:00:00 proplectic gregorian
units=times.units.split('since')[0]
basetime=pd.to_datetime(times.units.split('since')[-1])
dtimes = pd.to_timedelta(times[:],unit=str(units)[0])
tt=np.tile((basetime + dtimes),(H.shape[1],1)).T

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
ylim=(0,12)
xlim=None#(np.datetime64('2013-04-26T18:20:01'),np.datetime64('2013-04-26T20:59:00'))
xfmt = md.DateFormatter('%H:%M')
Ax = 2.0*(fr_rs_var['Attenuation_Hydrometeors'][:,0,:,0,0]+fr_rs_var['Attenuation_Atmosphere'][:,0,:,0])[:,::versus].cumsum(axis=1)[:,::versus]
Au = 2.0*(fr_rs_var['Attenuation_Hydrometeors'][:,0,:,1,0]+fr_rs_var['Attenuation_Atmosphere'][:,0,:,1])[:,::versus].cumsum(axis=1)[:,::versus]
Aa = 2.0*(fr_rs_var['Attenuation_Hydrometeors'][:,0,:,2,0]+fr_rs_var['Attenuation_Atmosphere'][:,0,:,2])[:,::versus].cumsum(axis=1)[:,::versus]
Aw = 2.0*(fr_rs_var['Attenuation_Hydrometeors'][:,0,:,3,0]+fr_rs_var['Attenuation_Atmosphere'][:,0,:,3])[:,::versus].cumsum(axis=1)[:,::versus]
Ag = 2.0*(fr_rs_var['Attenuation_Hydrometeors'][:,0,:,4,0]+fr_rs_var['Attenuation_Atmosphere'][:,0,:,4])[:,::versus].cumsum(axis=1)[:,::versus]
plot_variable(tt,H,Au,ax1,None,'height [km]','dB','Ku-band 2-way Attenuation',0,3,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Aa,ax2,None,'height [km]','dB','Ka-band 2-way Attenuation',0,3,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Aw,ax3,None,'height [km]','dB', 'W-band 2-way Attenuation',0,3,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Ag,ax4,'time','height [km]','dB', 'G-band 2-way Attenuation',0,10,xlim=xlim,ylim=ylim)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
xfmt = md.DateFormatter('%H:%M')
Zx = fr_rs_var['Ze'][:,0,:,0,0,0]
Zu = fr_rs_var['Ze'][:,0,:,1,0,0]
Za = fr_rs_var['Ze'][:,0,:,2,0,0]
Zw = fr_rs_var['Ze'][:,0,:,3,0,0]
Zg = fr_rs_var['Ze'][:,0,:,4,0,0]
plot_variable(tt,H,Zu,ax1,None,'height [km]','dBZ','Ku-band Ze',-20,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Za,ax2,None,'height [km]','dBZ','Ka-band Ze',-20,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Zw,ax3,None,'height [km]','dBZ', 'W-band Ze',-20,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Zg,ax4,'time','height [km]','dBZ', 'G-band Ze',-20,20,xlim=xlim,ylim=ylim)

ax4.xaxis.set_major_formatter(xfmt)
f.tight_layout(pad=0)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
DWRua = Zu-Za
DWRaw = Za-Zw
DWRag = Za-Zg
DWRwg = Zw-Zg
plot_variable(tt,H,DWRua,ax1,None,'height [km]','dB','DWR$_{Ku Ka}$',0,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,DWRaw,ax2,None,'height [km]','dB','DWR$_{Ka W}$',0,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,DWRag,ax3,None,'height [km]','dB','DWR$_{Ka G}$',0,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,DWRwg,ax4,'time','height [km]','dB','DWR$_{W G}$',0,20,xlim=xlim,ylim=ylim)
ax4.xaxis.set_major_formatter(xfmt)
f.tight_layout(pad=0)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
xfmt = md.DateFormatter('%H:%M')
Zx = fr_rs_var['Ze'][:,0,:,0,0,0]-Ax
Zu = fr_rs_var['Ze'][:,0,:,1,0,0]-Au
Za = fr_rs_var['Ze'][:,0,:,2,0,0]-Aa
Zw = fr_rs_var['Ze'][:,0,:,3,0,0]-Aw
Zg = fr_rs_var['Ze'][:,0,:,4,0,0]-Ag
plot_variable(tt,H,Zu,ax1,None,'height [km]','dBZ','Ku-band Z attenuated',-20,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Za,ax2,None,'height [km]','dBZ','Ka-band Z attenuated',-20,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Zw,ax3,None,'height [km]','dBZ', 'W-band Z attenuated',-20,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Zg,ax4,'time','height [km]','dBZ', 'G-band Z attenuated',-20,20,xlim=xlim,ylim=ylim)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
DWRua = Zu-Za
DWRaw = Za-Zw
DWRag = Za-Zg
DWRwg = Zw-Zg
plot_variable(tt,H,DWRua,ax1,None,'height [km]','dB','DWR$_{Ku Ka}$ attenuated',0,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,DWRaw,ax2,None,'height [km]','dB','DWR$_{Ka W}$ attenuated',0,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,DWRag,ax3,None,'height [km]','dB','DWR$_{Ka G}$ attenuated',0,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,DWRwg,ax4,'time','height [km]','dB','DWR$_{W G}$ attenuated',0,20,xlim=xlim,ylim=ylim)
ax4.xaxis.set_major_formatter(xfmt)
f.tight_layout(pad=0)


###############################################################################
# Plot hydrometeors

f,((ax1,ay1),(ax2,ay2),(ax3,ay3),(ax4,ay4),(ax6,ay6)) = plt.subplots(5,2,sharex=True,sharey=True)
var = fr_dt_var['QNI']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
var_dim_labels = var.dimensions
xvar = fr_dt_var[var_dim_labels[0]]
yvar = fr_dt_var[var_dim_labels[1]]
xval, yval = np.meshgrid(xvar[:], 0.001*yvar[:])
plot_variable(tt,H,vval,ax1,None,'height [Km]',var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = fr_dt_var['QNS']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
plot_variable(tt,H,vval,ax2,None,'height [Km]',var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = fr_dt_var['QNC']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
plot_variable(tt,H,vval,ax3,None,'height [Km]',var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = fr_dt_var['QNR']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
plot_variable(tt,H,vval,ax4,None,'height [Km]',var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = fr_dt_var['QNG']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
#plot_variable(tt,H,vval,ax5,None,'height [Km]',var.name+'  '+var.unit,var.long_name)
#var = fr_dt_var['QNH']
#vval = np.ma.masked_less(np.flip(var[:],1),0.1)
plot_variable(tt,H,vval,ax6,'time','height [Km]',var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)

var = fr_dt_var['QI']
vval = np.ma.masked_less(np.flip(var[:],1),1e-8)
plot_variable(tt,H,vval,ay1,None,None,var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = fr_dt_var['QS']
vval = np.ma.masked_less(np.flip(var[:],1),1e-8)
plot_variable(tt,H,vval,ay2,None,None,var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = fr_dt_var['QC']
vval = np.ma.masked_less(np.flip(var[:],1),1e-8)
plot_variable(tt,H,vval,ay3,None,None,var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = fr_dt_var['QR']
vval = np.ma.masked_less(np.flip(var[:],1),1e-8)
plot_variable(tt,H,vval,ay4,None,None,var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = fr_dt_var['QG']
vval = np.ma.masked_less(np.flip(var[:],1),1e-8)
#plot_variable(tt,H,vval,ay5,None,None,var.name+'  '+var.unit,var.long_name)
#var = fr_dt_var['QH']
#vval = np.ma.masked_less(np.flip(var[:],1),1e-7)
plot_variable(tt,H,vval,ay6,'time',None,var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
ax4.xaxis.set_major_formatter(xfmt)

diams = np.linspace(0.00005,0.0229,1000)
mass = lambda d: 0.038*d**2.0
vol = lambda d: (np.pi*d**3.0)/6.0
masses = mass(diams)
volumes = vol(diams)
Dm = lambda lam: ((masses*diams*diams*np.exp(-lam*diams)).sum()) / ((masses*diams*np.exp(-lam*diams)).sum())
D0 = lambda lam: ((volumes*diams*diams*np.exp(-lam*diams)).sum()) / ((volumes*diams*np.exp(-lam*diams)).sum())
vDm = np.vectorize(Dm)
vD0 = np.vectorize(D0)

lams = np.sqrt((0.228*np.flip(fr_dt_var['QNS'][:],1))/np.flip(fr_dt_var['QS'][:],1))
n0s = np.flip(fr_dt_var['QNS'][:],1)*lams**2.0
D0s = vD0(lams)
Dms = vDm(lams)
f,((ax1,ax2)) = plt.subplots(2,1,sharex=True,sharey=True)
plot_variable(tt,H,D0s,ax1,'time',None,'D0   [m]','Mean volume diameter',vmin=0.0,vmax=0.02,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Dms,ax2,'time',None,'Dm   [m]','Mean mass diameter',vmin=0.0,vmax=0.02,xlim=xlim,ylim=ylim)


###############################################################################
# Ny-Alesund case
#%%############################################################################
# Plot Z

H = ny_rs_var['height'][:,0,:]*0.001
times   = ny_dt_var['time'] # seconds since 2015-06-23 00:00:00 proplectic gregorian
units=times.units.split('since')[0]
basetime=pd.to_datetime(times.units.split('since')[-1])
dtimes = pd.to_timedelta(times[:],unit=str(units)[0])
tt=np.tile((basetime + dtimes),(H.shape[1],1)).T

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
ylim=(0,8)
xlim=(np.datetime64('2017-06-23T12:00:00'),np.datetime64('2017-06-23T23:00:00'))
xfmt = md.DateFormatter('%H:%M')
Ax = 2.0*(ny_rs_var['Attenuation_Hydrometeors'][:,0,:,0,0]+ny_rs_var['Attenuation_Atmosphere'][:,0,:,0])[:,::versus].cumsum(axis=1)[:,::versus]
Au = 2.0*(ny_rs_var['Attenuation_Hydrometeors'][:,0,:,1,0]+ny_rs_var['Attenuation_Atmosphere'][:,0,:,1])[:,::versus].cumsum(axis=1)[:,::versus]
Aa = 2.0*(ny_rs_var['Attenuation_Hydrometeors'][:,0,:,2,0]+ny_rs_var['Attenuation_Atmosphere'][:,0,:,2])[:,::versus].cumsum(axis=1)[:,::versus]
Aw = 2.0*(ny_rs_var['Attenuation_Hydrometeors'][:,0,:,3,0]+ny_rs_var['Attenuation_Atmosphere'][:,0,:,3])[:,::versus].cumsum(axis=1)[:,::versus]
Ag = 2.0*(ny_rs_var['Attenuation_Hydrometeors'][:,0,:,4,0]+ny_rs_var['Attenuation_Atmosphere'][:,0,:,4])[:,::versus].cumsum(axis=1)[:,::versus]
plot_variable(tt,H,Au,ax1,None,'height [km]','dB','Ku-band 2-way Attenuation',0,3,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Aa,ax2,None,'height [km]','dB','Ka-band 2-way Attenuation',0,3,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Aw,ax3,None,'height [km]','dB', 'W-band 2-way Attenuation',0,3,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Ag,ax4,'time','height [km]','dB', 'G-band 2-way Attenuation',0,10,xlim=xlim,ylim=ylim)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
xfmt = md.DateFormatter('%H:%M')
Zx = ny_rs_var['Ze'][:,0,:,0,0,0]
Zu = ny_rs_var['Ze'][:,0,:,1,0,0]
Za = ny_rs_var['Ze'][:,0,:,2,0,0]
Zw = ny_rs_var['Ze'][:,0,:,3,0,0]
Zg = ny_rs_var['Ze'][:,0,:,4,0,0]
plot_variable(tt,H,Zu,ax1,None,'height [km]','dBZ','Ku-band Ze',-20,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Za,ax2,None,'height [km]','dBZ','Ka-band Ze',-20,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Zw,ax3,None,'height [km]','dBZ', 'W-band Ze',-20,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Zg,ax4,'time','height [km]','dBZ', 'G-band Ze',-20,20,xlim=xlim,ylim=ylim)

ax4.xaxis.set_major_formatter(xfmt)
f.tight_layout(pad=0)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
DWRua = Zu-Za
DWRaw = Za-Zw
DWRag = Za-Zg
DWRwg = Zw-Zg
plot_variable(tt,H,DWRua,ax1,None,'height [km]','dB','DWR$_{Ku Ka}$',0,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,DWRaw,ax2,None,'height [km]','dB','DWR$_{Ka W}$',0,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,DWRag,ax3,None,'height [km]','dB','DWR$_{Ka G}$',0,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,DWRwg,ax4,'time','height [km]','dB','DWR$_{W G}$',0,20,xlim=xlim,ylim=ylim)
ax4.xaxis.set_major_formatter(xfmt)
f.tight_layout(pad=0)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
xfmt = md.DateFormatter('%H:%M')
Zx = ny_rs_var['Ze'][:,0,:,0,0,0]-Ax
Zu = ny_rs_var['Ze'][:,0,:,1,0,0]-Au
Za = ny_rs_var['Ze'][:,0,:,2,0,0]-Aa
Zw = ny_rs_var['Ze'][:,0,:,3,0,0]-Aw
Zg = ny_rs_var['Ze'][:,0,:,4,0,0]-Ag
plot_variable(tt,H,Zu,ax1,None,'height [km]','dBZ','Ku-band Z attenuated',-20,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Za,ax2,None,'height [km]','dBZ','Ka-band Z attenuated',-20,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Zw,ax3,None,'height [km]','dBZ', 'W-band Z attenuated',-20,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Zg,ax4,'time','height [km]','dBZ', 'G-band Z attenuated',-20,20,xlim=xlim,ylim=ylim)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
DWRua = Zu-Za
DWRaw = Za-Zw
DWRag = Za-Zg
DWRwg = Zw-Zg
plot_variable(tt,H,DWRua,ax1,None,'height [km]','dB','DWR$_{Ku Ka}$ attenuated',0,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,DWRaw,ax2,None,'height [km]','dB','DWR$_{Ka W}$ attenuated',0,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,DWRag,ax3,None,'height [km]','dB','DWR$_{Ka G}$ attenuated',0,20,xlim=xlim,ylim=ylim)
plot_variable(tt,H,DWRwg,ax4,'time','height [km]','dB','DWR$_{W G}$ attenuated',0,20,xlim=xlim,ylim=ylim)
ax4.xaxis.set_major_formatter(xfmt)
f.tight_layout(pad=0)


###############################################################################
# Plot hydrometeors

f,((ax1,ay1),(ax2,ay2),(ax3,ay3),(ax4,ay4),(ax6,ay6)) = plt.subplots(5,2,sharex=True,sharey=True)
var = ny_dt_var['QNI']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
var_dim_labels = var.dimensions
xvar = ny_dt_var[var_dim_labels[0]]
yvar = ny_dt_var[var_dim_labels[1]]
xval, yval = np.meshgrid(xvar[:], 0.001*yvar[:])
plot_variable(tt,H,vval,ax1,None,'height [Km]',var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = ny_dt_var['QNS']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
plot_variable(tt,H,vval,ax2,None,'height [Km]',var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = ny_dt_var['QNC']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
plot_variable(tt,H,vval,ax3,None,'height [Km]',var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = ny_dt_var['QNR']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
plot_variable(tt,H,vval,ax4,None,'height [Km]',var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = ny_dt_var['QNG']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
#plot_variable(tt,H,vval,ax5,None,'height [Km]',var.name+'  '+var.unit,var.long_name)
#var = ny_dt_var['QNH']
#vval = np.ma.masked_less(np.flip(var[:],1),0.1)
plot_variable(tt,H,vval,ax6,'time','height [Km]',var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)

var = ny_dt_var['QI']
vval = np.ma.masked_less(np.flip(var[:],1),1e-8)
plot_variable(tt,H,vval,ay1,None,None,var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = ny_dt_var['QS']
vval = np.ma.masked_less(np.flip(var[:],1),1e-8)
plot_variable(tt,H,vval,ay2,None,None,var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = ny_dt_var['QC']
vval = np.ma.masked_less(np.flip(var[:],1),1e-8)
plot_variable(tt,H,vval,ay3,None,None,var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = ny_dt_var['QR']
vval = np.ma.masked_less(np.flip(var[:],1),1e-8)
plot_variable(tt,H,vval,ay4,None,None,var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
var = ny_dt_var['QG']
vval = np.ma.masked_less(np.flip(var[:],1),1e-8)
#plot_variable(tt,H,vval,ay5,None,None,var.name+'  '+var.unit,var.long_name)
#var = ny_dt_var['QH']
#vval = np.ma.masked_less(np.flip(var[:],1),1e-7)
plot_variable(tt,H,vval,ay6,'time',None,var.name+'  '+var.unit,var.long_name,xlim=xlim,ylim=ylim)
ax4.xaxis.set_major_formatter(xfmt)

lams = np.sqrt((0.228*np.flip(ny_dt_var['QNS'][:],1))/np.flip(ny_dt_var['QS'][:],1))
n0s = np.flip(ny_dt_var['QNS'][:],1)*lams**2.0
D0s = vD0(lams)
Dms = vDm(lams)
f,((ax1,ax2)) = plt.subplots(2,1,sharex=True,sharey=True)
plot_variable(tt,H,D0s,ax1,'time',None,'D0   [m]','Mean volume diameter',vmin=0.0,vmax=0.02,xlim=xlim,ylim=ylim)
plot_variable(tt,H,Dms,ax2,'time',None,'Dm   [m]','Mean mass diameter',vmin=0.0,vmax=0.02,xlim=xlim,ylim=ylim)

#%%############################################################################
# 24 November case
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
Ax = 2.0*(ml_rs_var['Attenuation_Hydrometeors'][:,0,:,0,0]+ml_rs_var['Attenuation_Atmosphere'][:,0,:,0])[:,::versus].cumsum(axis=1)[:,::versus]
Au = 2.0*(ml_rs_var['Attenuation_Hydrometeors'][:,0,:,1,0]+ml_rs_var['Attenuation_Atmosphere'][:,0,:,1])[:,::versus].cumsum(axis=1)[:,::versus]
Aa = 2.0*(ml_rs_var['Attenuation_Hydrometeors'][:,0,:,2,0]+ml_rs_var['Attenuation_Atmosphere'][:,0,:,2])[:,::versus].cumsum(axis=1)[:,::versus]
Aw = 2.0*(ml_rs_var['Attenuation_Hydrometeors'][:,0,:,3,0]+ml_rs_var['Attenuation_Atmosphere'][:,0,:,3])[:,::versus].cumsum(axis=1)[:,::versus]
Ag = 2.0*(ml_rs_var['Attenuation_Hydrometeors'][:,0,:,4,0]+ml_rs_var['Attenuation_Atmosphere'][:,0,:,4])[:,::versus].cumsum(axis=1)[:,::versus]
plot_variable(tt,H,Au,ax1,None,'height [km]','dB','Ku-band 2-way Attenuation',0,3,ylim=ylim)
plot_variable(tt,H,Aa,ax2,None,'height [km]','dB','Ka-band 2-way Attenuation',0,3,ylim=ylim)
plot_variable(tt,H,Aw,ax3,None,'height [km]','dB', 'W-band 2-way Attenuation',0,3,ylim=ylim)
plot_variable(tt,H,Ag,ax4,'time','height [km]','dB', 'G-band 2-way Attenuation',0,10,ylim=ylim)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
ylim=(0,15)
xfmt = md.DateFormatter('%H:%M')
Zx = ml_rs_var['Ze'][:,0,:,0,0,0]
Zu = ml_rs_var['Ze'][:,0,:,1,0,0]
Za = ml_rs_var['Ze'][:,0,:,2,0,0]
Zw = ml_rs_var['Ze'][:,0,:,3,0,0]
Zg = ml_rs_var['Ze'][:,0,:,4,0,0]
plot_variable(tt,H,Zu,ax1,None,'height [km]','dBZ','Ku-band Ze',-20,20,ylim=ylim)
plot_variable(tt,H,Za,ax2,None,'height [km]','dBZ','Ka-band Ze',-20,20,ylim=ylim)
plot_variable(tt,H,Zw,ax3,None,'height [km]','dBZ', 'W-band Ze',-20,20,ylim=ylim)
plot_variable(tt,H,Zg,ax4,'time','height [km]','dBZ', 'G-band Ze',-20,20,ylim=ylim)

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

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
ylim=(0,15)
xfmt = md.DateFormatter('%H:%M')
Zx = ml_rs_var['Ze'][:,0,:,0,0,0]-Ax
Zu = ml_rs_var['Ze'][:,0,:,1,0,0]-Au
Za = ml_rs_var['Ze'][:,0,:,2,0,0]-Aa
Zw = ml_rs_var['Ze'][:,0,:,3,0,0]-Aw
Zg = ml_rs_var['Ze'][:,0,:,4,0,0]-Ag
plot_variable(tt,H,Zu,ax1,None,'height [km]','dBZ','Ku-band Z attenuated',-20,20,ylim=ylim)
plot_variable(tt,H,Za,ax2,None,'height [km]','dBZ','Ka-band Z attenuated',-20,20,ylim=ylim)
plot_variable(tt,H,Zw,ax3,None,'height [km]','dBZ', 'W-band Z attenuated',-20,20,ylim=ylim)
plot_variable(tt,H,Zg,ax4,'time','height [km]','dBZ', 'G-band Z attenuated',-20,20,ylim=ylim)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
DWRua = Zu-Za
DWRaw = Za-Zw
DWRag = Za-Zg
DWRwg = Zw-Zg
plot_variable(tt,H,DWRua,ax1,None,'height [km]','dB','DWR$_{Ku Ka}$ attenuated',0,20,ylim=ylim)
plot_variable(tt,H,DWRaw,ax2,None,'height [km]','dB','DWR$_{Ka W}$ attenuated',0,20,ylim=ylim)
plot_variable(tt,H,DWRag,ax3,None,'height [km]','dB','DWR$_{Ka G}$ attenuated',0,20,ylim=ylim)
plot_variable(tt,H,DWRwg,ax4,'time','height [km]','dB','DWR$_{W G}$ attenuated',0,20,ylim=ylim)
ax4.xaxis.set_major_formatter(xfmt)
f.tight_layout(pad=0)


###############################################################################
# Plot hydrometeors

f,((ax1,ay1),(ax2,ay2),(ax3,ay3),(ax4,ay4),(ax6,ay6)) = plt.subplots(5,2,sharex=True,sharey=True)
var = ml_dt_var['QNI']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
var_dim_labels = var.dimensions
xvar = ml_dt_var[var_dim_labels[0]]
yvar = ml_dt_var[var_dim_labels[1]]
xval, yval = np.meshgrid(xvar[:], 0.001*yvar[:])
plot_variable(tt,H,vval,ax1,None,'height [Km]',var.name+'  '+var.unit,var.long_name,ylim=ylim)
var = ml_dt_var['QNS']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
plot_variable(tt,H,vval,ax2,None,'height [Km]',var.name+'  '+var.unit,var.long_name,ylim=ylim)
var = ml_dt_var['QNC']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
plot_variable(tt,H,vval,ax3,None,'height [Km]',var.name+'  '+var.unit,var.long_name,ylim=ylim)
var = ml_dt_var['QNR']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
plot_variable(tt,H,vval,ax4,None,'height [Km]',var.name+'  '+var.unit,var.long_name,ylim=ylim)
var = ml_dt_var['QNG']
vval = np.ma.masked_less(np.flip(var[:],1),0.1)
#plot_variable(tt,H,vval,ax5,None,'height [Km]',var.name+'  '+var.unit,var.long_name,ylim=ylim)
#var = ml_dt_var['QNH']
#vval = np.ma.masked_less(np.flip(var[:],1),0.1)
plot_variable(tt,H,vval,ax6,'time','height [Km]',var.name+'  '+var.unit,var.long_name,ylim=ylim)

var = ml_dt_var['QI']
vval = np.ma.masked_less(np.flip(var[:],1),1e-7)
plot_variable(tt,H,vval,ay1,None,None,var.name+'  '+var.unit,var.long_name,ylim=ylim)
var = ml_dt_var['QS']
vval = np.ma.masked_less(np.flip(var[:],1),1e-7)
plot_variable(tt,H,vval,ay2,None,None,var.name+'  '+var.unit,var.long_name,ylim=ylim)
var = ml_dt_var['QC']
vval = np.ma.masked_less(np.flip(var[:],1),1e-7)
plot_variable(tt,H,vval,ay3,None,None,var.name+'  '+var.unit,var.long_name,ylim=ylim)
var = ml_dt_var['QR']
vval = np.ma.masked_less(np.flip(var[:],1),1e-7)
plot_variable(tt,H,vval,ay4,None,None,var.name+'  '+var.unit,var.long_name,ylim=ylim)
var = ml_dt_var['QG']
vval = np.ma.masked_less(np.flip(var[:],1),1e-7)
#plot_variable(tt,H,vval,ay5,None,None,var.name+'  '+var.unit,var.long_name,ylim=ylim)
#var = ml_dt_var['QH']
#vval = np.ma.masked_less(np.flip(var[:],1),1e-7)
plot_variable(tt,H,vval,ay6,'time',None,var.name+'  '+var.unit,var.long_name,ylim=ylim)
ax4.xaxis.set_major_formatter(xfmt)

lams = np.sqrt((0.228*np.flip(ml_dt_var['QNS'][:],1))/np.flip(ml_dt_var['QS'][:],1))
n0s = np.flip(ml_dt_var['QNS'][:],1)*lams**2.0
D0s = vD0(lams)
Dms = vDm(lams)
f,((ax1,ax2)) = plt.subplots(2,1,sharex=True,sharey=True)
plot_variable(tt,H,D0s,ax1,'time',None,'D0   [m]','Mean volume diameter',vmin=0.0,vmax=0.02,ylim=ylim)
plot_variable(tt,H,Dms,ax2,'time',None,'Dm   [m]','Mean mass diameter',vmin=0.0,vmax=0.02,ylim=ylim)


###############################################################################
#%%#      MIXED PHASE CASE
###############################################################################
# Plot hydrometeors

var = mp_dt_var['QNI']
vval = np.ma.masked_less(var[:],0.01) #var[:].T 
var_dim_labels = var.dimensions
xvar = mp_dt_var[var_dim_labels[1]]
yvar = mp_dt_var['height_m']
xval, yval = np.meshgrid(xvar[:], 0.001*yvar[:])

f,((ax1,ax2),(ax3,ax4),(ax5,ax6),(ax7,ax8)) = plt.subplots(4,2,sharex=True,sharey=True)
ylim = (0,5)
plot_variable(xval,yval,vval,ax1,None,'height [Km]','QNI kg-1','number concentration ice',ylim=ylim)
var = mp_dt_var['QNS']
vval = np.ma.masked_less(var[:],0.01) #var[:].T
plot_variable(xval,yval,vval,ax3,None,'height [Km]','QNS kg-1','number concentration snow',ylim=ylim)
var = mp_dt_var['QNC']
vval = np.ma.masked_less(var[:],0.01) #var[:].T
plot_variable(xval,yval,vval,ax5,None,'height [Km]','QNC kg-1','number concentration cloud droplets',ylim=ylim)
var = mp_dt_var['QNG']
vval = np.ma.masked_less(var[:],0.01) #var[:].T
plot_variable(xval,yval,vval,ax7,xvar.name,'height [Km]','QNG kg-1','number concentration graupel',ylim=ylim)

var = mp_dt_var['QI']
vval = np.ma.masked_less(var[:],1e-8) #var[:].T
plot_variable(xval,yval,vval,ax2,None,None,'QI  kg kg-1','specific ice content' ,ylim=ylim)
var = mp_dt_var['QS']
vval = np.ma.masked_less(var[:],1e-8) #var[:].T
plot_variable(xval,yval,vval,ax4,None,None,'QS  kg kg-1','snow mixing ratio' ,ylim=ylim)
var = mp_dt_var['QC']
vval = np.ma.masked_less(var[:],1e-8) #var[:].T
plot_variable(xval,yval,vval,ax6,None,None,'QC  kg kg-1','cloud droplets mixing ratio',ylim=ylim)
var = mp_dt_var['QG']
vval = np.ma.masked_less(var[:],1e-8) #var[:].T
plot_variable(xval,yval,vval,ax8,xvar.name,None,'QG  kg kg-1','graupel mixing ratio',ylim=ylim)
f.tight_layout(pad=0)

lams = np.sqrt((0.228*mp_dt_var['QNS'][:])/(mp_dt_var['QS'][:]))
n0s = np.flip(mp_dt_var['QNS'][:],1)*lams**2.0
D0s = vD0(lams)
Dms = vDm(lams)
f,((ax1,ax2)) = plt.subplots(2,1,sharex=True,sharey=True)
plot_variable(xval,yval,D0s,ax1,xvar.name,None,'D0   [m]','Mean volume diameter',vmin=0.0,vmax=0.02,ylim=ylim)
plot_variable(xval,yval,Dms,ax2,xvar.name,None,'Dm   [m]','Mean mass diameter',vmin=0.0,vmax=0.02,ylim=ylim)

###############################################################################
# Plot Z
f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
ylim=(0,5)
xfmt = md.DateFormatter('%H:%M')
Ax = 2.0*(mp_rs_var['Attenuation_Hydrometeors'][:,0,:,0,0]+mp_rs_var['Attenuation_Atmosphere'][:,0,:,0])[:,::versus].cumsum(axis=1)[:,::versus]
Au = 2.0*(mp_rs_var['Attenuation_Hydrometeors'][:,0,:,1,0]+mp_rs_var['Attenuation_Atmosphere'][:,0,:,1])[:,::versus].cumsum(axis=1)[:,::versus]
Aa = 2.0*(mp_rs_var['Attenuation_Hydrometeors'][:,0,:,2,0]+mp_rs_var['Attenuation_Atmosphere'][:,0,:,2])[:,::versus].cumsum(axis=1)[:,::versus]
Aw = 2.0*(mp_rs_var['Attenuation_Hydrometeors'][:,0,:,3,0]+mp_rs_var['Attenuation_Atmosphere'][:,0,:,3])[:,::versus].cumsum(axis=1)[:,::versus]
Ag = 2.0*(mp_rs_var['Attenuation_Hydrometeors'][:,0,:,4,0]+mp_rs_var['Attenuation_Atmosphere'][:,0,:,4])[:,::versus].cumsum(axis=1)[:,::versus]
plot_variable(xval,yval,np.flip(Au.T,0),ax1,None,'height [km]','dB','Ku-band 2-way Attenuation',0,5,ylim=ylim)
plot_variable(xval,yval,np.flip(Aa.T,0),ax2,None,'height [km]','dB','Ka-band 2-way Attenuation',0,5,ylim=ylim)
plot_variable(xval,yval,np.flip(Aw.T,0),ax3,None,'height [km]','dB', 'W-band 2-way Attenuation',0,20,ylim=ylim)
plot_variable(xval,yval,np.flip(Ag.T,0),ax4,'lat','height [km]','dB', 'G-band 2-way Attenuation',0,20,ylim=ylim)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
Zx = mp_rs_var['Ze'][:,0,:,0,0,0]
Zu = mp_rs_var['Ze'][:,0,:,1,0,0]
Za = mp_rs_var['Ze'][:,0,:,2,0,0]
Zw = mp_rs_var['Ze'][:,0,:,3,0,0]
Zg = mp_rs_var['Ze'][:,0,:,4,0,0]
plot_variable(xval,yval,np.flip(Zu.T,0),ax1,None,'height [km]','dBZ','Ku-band Ze',-20,20,ylim=ylim)
plot_variable(xval,yval,np.flip(Za.T,0),ax2,None,'height [km]','dBZ','Ka-band Ze',-20,20,ylim=ylim)
plot_variable(xval,yval,np.flip(Zw.T,0),ax3,None,'height [km]','dBZ', 'W-band Ze',-20,20,ylim=ylim)
plot_variable(xval,yval,np.flip(Zg.T,0),ax4,'lat','height [km]','dBZ', 'G-band Ze',-20,20,ylim=ylim)
f.tight_layout(pad=0)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
DWRua = Zu-Za
DWRaw = Za-Zw
DWRag = Za-Zg
DWRwg = Zw-Zg
plot_variable(xval,yval,np.flip(DWRua.T,0),ax1,None,'height [km]','dB','DWR$_{Ku Ka}$',0,20,ylim=ylim)
plot_variable(xval,yval,np.flip(DWRaw.T,0),ax2,None,'height [km]','dB','DWR$_{Ka W}$',0,20,ylim=ylim)
plot_variable(xval,yval,np.flip(DWRag.T,0),ax3,None,'height [km]','dB', 'DWR$_{Ka G}$',0,20,ylim=ylim)
plot_variable(xval,yval,np.flip(DWRwg.T,0),ax4,'lat','height [km]','dB', 'DWR$_{W G}$',0,20,ylim=ylim)
f.tight_layout(pad=0)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
Zx = mp_rs_var['Ze'][:,0,:,0,0,0]-Ax
Zu = mp_rs_var['Ze'][:,0,:,1,0,0]-Au
Za = mp_rs_var['Ze'][:,0,:,2,0,0]-Aa
Zw = mp_rs_var['Ze'][:,0,:,3,0,0]-Aw
Zg = mp_rs_var['Ze'][:,0,:,4,0,0]-Ag
plot_variable(xval,yval,np.flip(Zu.T,0),ax1,None,'height [km]','dBZ','Ku-band Ze attenuated',-20,20,ylim=ylim)
plot_variable(xval,yval,np.flip(Za.T,0),ax2,None,'height [km]','dBZ','Ka-band Ze attenuated',-20,20,ylim=ylim)
plot_variable(xval,yval,np.flip(Zw.T,0),ax3,None,'height [km]','dBZ', 'W-band Ze attenuated',-20,20,ylim=ylim)
plot_variable(xval,yval,np.flip(Zg.T,0),ax4,'lat','height [km]','dBZ', 'G-band Ze attenuated',-20,20,ylim=ylim)
f.tight_layout(pad=0)

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=True)
DWRua = Zu-Za
DWRaw = Za-Zw
DWRag = Za-Zg
DWRwg = Zw-Zg
plot_variable(xval,yval,np.flip(DWRua.T,0),ax1,None,'height [km]','dB','DWR$_{Ku Ka} attenuated$',0,20,ylim=ylim)
plot_variable(xval,yval,np.flip(DWRaw.T,0),ax2,None,'height [km]','dB','DWR$_{Ka W}$ attenuated',0,20,ylim=ylim)
plot_variable(xval,yval,np.flip(DWRag.T,0),ax3,None,'height [km]','dB', 'DWR$_{Ka G}$ attenuated',0,20,ylim=ylim)
plot_variable(xval,yval,np.flip(DWRwg.T,0),ax4,'lat','height [km]','dB', 'DWR$_{W G}$ attenuated',0,20,ylim=ylim)
f.tight_layout(pad=0)

