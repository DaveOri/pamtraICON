import netCDF4
import matplotlib.pyplot as plt
plt.close('all')
import matplotlib.dates as md
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
from glob import glob
import os

patch = '001'

iconfile = '/net/junta/schemann/ac3prop/dendrite/METEOGRAM_patch'+patch+'_awipev.nc'
datafile = '/data/optimice/pamtra_runs/nyalesund/20170602/METEOGRAM_02062017_awipev_patch'+patch+'_dendrite.nc'
plotpath = '/data/optimice/pamtra_runs/nyalesund/20170602/20170602dendrite_P'+patch+'_'

iconfile = '/net/junta/schemann/ac3prop/default/METEOGRAM_patch'+patch+'_awipev.nc'
datafile = '/data/optimice/pamtra_runs/nyalesund/20170602/METEOGRAM_02062017_awipev_patch'+patch+'_icecosmo.nc'
plotpath = '/data/optimice/pamtra_runs/nyalesund/20170602/20170602icecosmo_P'+patch+'_'

#iconfile = '/net/junta/schemann/ac3prop/icehex/METEOGRAM_patch'+patch+'_awipev.nc'
#datafile = '/data/optimice/pamtra_runs/nyalesund/20170602/METEOGRAM_02062017_awipev_patch'+patch+'_icehex.nc'
#plotpath = '/data/optimice/pamtra_runs/nyalesund/20170602/20170602icehex_P'+patch+'_'

figsize21 = (18,12)
figsize31 = (18,18)
figsize41 = (18,24)
figsize51 = (18,30)

figsize21 = (8,6.4)
figsize31 = (8,9.6)
figsize41 = (8,12.8)
figsize51 = (8,16)

versus = -1 # Top Down
versus =  1 # Bottom Up
#xfmt = md.DateFormatter('%m-%d %H')
xfmt = md.DateFormatter('%H')
#ylim=(0,3000)
ylim=(0,1500)
xDataLim = -2

dpi = 150

def plot_variable(x,y,v,axes,
                  xlab=None,ylab=None,vlab=None,title=None,
                  vmin=None,vmax=None,xlim=None,ylim=None,
                  cmap='jet', **kwargs):
    mesh = axes.pcolormesh(x,y,v,vmin=vmin,vmax=vmax,cmap=cmap, **kwargs)
    if title is not None:
        axes.text(0.1,0.9,title,transform=axes.transAxes,
                  weight='black', size='x-large',
                  bbox=dict(facecolor='white'))
    plt.colorbar(mesh,label=vlab,ax=axes)
    if xlab is not None:
        axes.set_xlabel(xlab)
    if ylab is not None:
        axes.set_ylabel(ylab)
    axes.set_xlim(xlim)
    axes.set_ylim(ylim)

runDataset = netCDF4.Dataset(datafile)
runVars = runDataset.variables

iconvars = netCDF4.Dataset(iconfile).variables
Hicon = iconvars['height_2'][:]
Ticon = netCDF4.num2date(iconvars['time'][:], iconvars['time'].units)
qc = iconvars['QC'][:].T
qg = iconvars['QG'][:].T
qs = iconvars['QS'][:].T

f,((ax1,ax2,ax3)) = plt.subplots(3, 1, sharex=False, figsize=figsize31)
plot_variable(Ticon,Hicon,qc,ax1,None,'height [m]','[kg/kg]','Cloud droplets mixing ratio',ylim=ylim)
plot_variable(Ticon,Hicon,qg,ax2,None,'height [m]','[kg/kg]','Graupel mixing ratio',ylim=ylim)
plot_variable(Ticon,Hicon,qs,ax3,'time  UTC','height [m]','kg/kg', 'Snow mixing ratio',ylim=ylim)
ax1.xaxis.set_major_formatter(xfmt)
ax2.xaxis.set_major_formatter(xfmt)
ax3.xaxis.set_major_formatter(xfmt)
ax1.grid()
ax2.grid()
ax3.grid()
ax1.xaxis.set_major_formatter(xfmt)
ax2.xaxis.set_major_formatter(xfmt)
ax3.xaxis.set_major_formatter(xfmt)
f.tight_layout(pad=0)
f.savefig(plotpath+'mixing_ratios'+'.png', dpi=dpi, bbox_inches='tight')

# Extract Heights
H = (runVars['height'][:,0,:])[:xDataLim,:]
# Reshape times
ttt = pd.to_datetime(runVars['datatime'][:,0],unit='s')
tt = (np.tile(ttt,(H.shape[1],1)).T)[:xDataLim,:]
print(tt.shape, H.shape)
# Extract attenuation
#ax = 2.0*(runVars['Attenuation_Hydrometeors'][:,0,:,0,0] + 
#    runVars['Attenuation_Atmosphere'][:,0,:,0])
#Ax = ax[:,::versus].cumsum(axis=1)[:,::versus][:xDataLim,:]
#au = 2.0*(runVars['Attenuation_Hydrometeors'][:,0,:,1,0] + 
#    runVars['Attenuation_Atmosphere'][:,0,:,1])
#Au = au[:,::versus].cumsum(axis=1)[:,::versus][:xDataLim,:]
aa = 2.0*(runVars['Attenuation_Hydrometeors'][:,0,:,0,0] + 
    runVars['Attenuation_Atmosphere'][:,0,:,0])
Aa = aa[:,::versus].cumsum(axis=1)[:,::versus][:xDataLim,:]
aw = 2.0*(runVars['Attenuation_Hydrometeors'][:,0,:,1,0] + 
    runVars['Attenuation_Atmosphere'][:,0,:,1])
Aw = aw[:,::versus].cumsum(axis=1)[:,::versus][:xDataLim,:]
#ag = 2.0*(runVars['Attenuation_Hydrometeors'][:,0,:,4,0] + 
#    runVars['Attenuation_Atmosphere'][:,0,:,4])
#Ag = ag[:,::versus].cumsum(axis=1)[:,::versus][:xDataLim,:]

# Plot Attenuation
#f,((ax1,ax2,ax3,ax4)) = plt.subplots(4, 1, sharex=False, figsize=figsize41)
f,((ax1,ax2)) = plt.subplots(2, 1, sharex=False, figsize=figsize21)
#plot_variable(tt,H,Ax,ax1,None,'height [m]','dB','X-band 2-way Attenuation',0,1,ylim=ylim)
plot_variable(tt,H,Aa,ax1,None,'height [m]','dB','Ka-band 2-way Attenuation',0,5,ylim=ylim)
#plot_variable(tt,H,Aw,ax3,None,'height [m]','dB', 'W-band 2-way Attenuation',0,15,ylim=ylim)
plot_variable(tt,H,Aw,ax2,'time  UTC','height [m]','dB', 'W-band 2-way Attenuation',0,15,ylim=ylim)
#plot_variable(tt,H,Ag,ax4,'time  UTC','height [m]','dB', 'G-band 2-way Attenuation',0,15,ylim=ylim)
#ax1.set_title('Ka-band')
#ax2.set_title('W-band')
##ax3.set_title('W-band')
##ax4.set_title('G-band')
ax1.xaxis.set_major_formatter(xfmt)
ax2.xaxis.set_major_formatter(xfmt)
#ax3.xaxis.set_major_formatter(xfmt)
#ax4.xaxis.set_major_formatter(xfmt)
ax1.grid()
ax2.grid()
#ax3.grid()
#ax4.grid()
f.tight_layout(pad=0)
f.savefig(plotpath+'attenuation'+'.png', dpi=dpi, bbox_inches='tight')

# Extract Ze
#Zex = runVars['Ze'][:,0,:,0,0,0][:xDataLim,:]
#Zeu = runVars['Ze'][:,0,:,1,0,0][:xDataLim,:]
Zea = runVars['Ze'][:,0,:,0,0,0][:xDataLim,:]
Zew = runVars['Ze'][:,0,:,1,0,0][:xDataLim,:]
#Zeg = runVars['Ze'][:,0,:,4,0,0][:xDataLim,:]

# make DWRs and plot
#DWRxa = Zex-Zea
DWRaw = Zea-Zew
#DWRag = Zea-Zeg
#DWRwg = Zew-Zeg

# Plot Ze
#f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=False,figsize=figsize41)
f,((ax1,ax2,ax3)) = plt.subplots(3,1,sharex=False,figsize=figsize31)
#plot_variable(tt,H,Zex,ax1,None,'height [m]','dBZ','X-band Ze',-50,-10,ylim=ylim)
plot_variable(tt,H,Zea,ax1,None,'height [m]','dBZ', 'Ka-band Ze',-50,-10,ylim=ylim)
plot_variable(tt,H,Zew,ax2,None,'height [m]','dBZ', 'W-band Ze',-50,-10,ylim=ylim)
plot_variable(tt,H,DWRaw,ax3,'time  UTC','height [m]','dB','DWR$_{Ka W}$',-5,20, ylim=ylim,cmap='nipy_spectral')
#plot_variable(tt,H,Zeg,ax4,'time  UTC','height [m]','dBZ', 'G-band Ze',-50,-10,ylim=ylim)
ax1.xaxis.set_major_formatter(xfmt)
ax2.xaxis.set_major_formatter(xfmt)
ax3.xaxis.set_major_formatter(xfmt)
#ax4.xaxis.set_major_formatter(xfmt)
#ax1.set_title('Ka-band')
#ax2.set_title('W-band')
#ax3.set_title('DWRkw-band')
##ax4.set_title('G-band')
ax1.grid()
ax2.grid()
ax3.grid()
#ax4.grid()
f.tight_layout(pad=0)
f.savefig(plotpath+'Ze_DWR'+'.png',dpi=dpi, bbox_inches='tight')



# f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,1,sharex=False,figsize=figsize41)
# plot_variable(tt,H,DWRxa,ax1,None,'height [m]','dB','DWR$_{X Ka}$',-5,20, ylim=ylim,cmap='nipy_spectral')
# plot_variable(tt,H,DWRaw,ax2,None,'height [m]','dB','DWR$_{Ka W}$',-5,20, ylim=ylim,cmap='nipy_spectral')
# plot_variable(tt,H,DWRag,ax3,None,'height [m]','dB','DWR$_{Ka G}$',-5,20, ylim=ylim,cmap='nipy_spectral')
# plot_variable(tt,H,DWRwg,ax4,'time  UTC','height [m]','dB','DWR$_{W G}$',-5,10, ylim=ylim,cmap='nipy_spectral')
# ax1.xaxis.set_major_formatter(xfmt)
# ax2.xaxis.set_major_formatter(xfmt)
# ax3.xaxis.set_major_formatter(xfmt)
# ax4.xaxis.set_major_formatter(xfmt)
# #ax1.set_title('X-Ka')
# #ax2.set_title('Ka-W')
# #ax3.set_title('Ka-G')
# #ax4.set_title('W-G')
# ax1.grid()
# ax2.grid()
# ax3.grid()
# ax4.grid()
# f.tight_layout(pad=0)
# f.savefig(plotpath+'DWRe'+'.png', dpi=dpi, bbox_inches='tight')

# make attenuated Z and DWRs and respective plots
#Zx = Zex-Ax
#Zu = Zeu-Au
Za = Zea-Aa
Zw = Zew-Aw
#Zg = Zeg-Ag

#DWRxa = Zx-Za
DWRaw = Za-Zw
#DWRag = Za-Zg
#DWRwg = Zw-Zg
f,((ax1,ax2,ax3)) = plt.subplots(3,1,sharex=False,figsize=figsize31)
#plot_variable(tt,H,Zex,ax1,None,'height [m]','dBZ','X-band Ze',-50,-10,ylim=ylim)
plot_variable(tt,H,Za,ax1,None,'height [m]','dBZ', 'Ka-band Ze',-50,-10,ylim=ylim)
plot_variable(tt,H,Zw,ax2,None,'height [m]','dBZ', 'W-band Ze',-50,-10,ylim=ylim)
plot_variable(tt,H,DWRaw,ax3,'time  UTC','height [m]','dB','DWR$_{Ka W}$',-5,20, ylim=ylim,cmap='nipy_spectral')
#plot_variable(tt,H,Zeg,ax4,'time  UTC','height [m]','dBZ', 'G-band Ze',-50,-10,ylim=ylim)
ax1.xaxis.set_major_formatter(xfmt)
ax2.xaxis.set_major_formatter(xfmt)
ax3.xaxis.set_major_formatter(xfmt)
#ax4.xaxis.set_major_formatter(xfmt)
#ax1.set_title('Ka-band attenuated')
#ax2.set_title('W-band attenuated')
#ax3.set_title('DWRkw-band attenuated')
##ax4.set_title('G-band')
ax1.grid()
ax2.grid()
ax3.grid()
#ax4.grid()
f.tight_layout(pad=0)
f.savefig(plotpath+'Ze_DWR_attenuated'+'.png',dpi=dpi, bbox_inches='tight')

# Extract MDV
#MDVx = -runVars['Radar_MeanDopplerVel'][:,0,:,0,0,0][:xDataLim,:]
#MDVu = -runVars['Radar_MeanDopplerVel'][:,0,:,1,0,0][:xDataLim,:]
MDVa = -runVars['Radar_MeanDopplerVel'][:,0,:,0,0,0][:xDataLim,:]
MDVw = -runVars['Radar_MeanDopplerVel'][:,0,:,1,0,0][:xDataLim,:]
#MDVg = -runVars['Radar_MeanDopplerVel'][:,0,:,4,0,0][:xDataLim,:]

# Extract SW
#SWx = runVars['Radar_SpectrumWidth'][:,0,:,0,0,0][:xDataLim,:]
#SWu = runVars['Radar_SpectrumWidth'][:,0,:,1,0,0][:xDataLim,:]
SWa = runVars['Radar_SpectrumWidth'][:,0,:,0,0,0][:xDataLim,:]
SWw = runVars['Radar_SpectrumWidth'][:,0,:,1,0,0][:xDataLim,:]
#SWg = runVars['Radar_SpectrumWidth'][:,0,:,4,0,0][:xDataLim,:]

# Plot MDV
f,((ax1,ax2)) = plt.subplots(2, 1, sharex=False, figsize=figsize21)
plot_variable(tt,H,MDVa,ax1,None,'height [m]','m/s','Ka-band Mean Doppler Velocity',0,2,ylim=ylim)
plot_variable(tt,H,MDVw,ax2,'time  UTC','height [m]','m/s', 'W-band Mean Doppler Velocity',0,2,ylim=ylim)
#ax1.set_title('Ka-band')
#ax2.set_title('W-band')
ax1.xaxis.set_major_formatter(xfmt)
ax2.xaxis.set_major_formatter(xfmt)
ax1.grid()
ax2.grid()
f.tight_layout(pad=0)
f.savefig(plotpath+'MDV'+'.png', dpi=dpi, bbox_inches='tight')

# Plot SW
f,((ax1,ax2)) = plt.subplots(2, 1, sharex=False, figsize=figsize21)
plot_variable(tt,H,SWa,ax1,None,'height [m]','m/s','Ka-band Spectrum Width',0,2,ylim=ylim,cmap='nipy_spectral')
plot_variable(tt,H,SWw,ax2,'time  UTC','height [m]','m/s', 'W-band Spectrum Width',0,2,ylim=ylim,cmap='nipy_spectral')
#ax1.set_title('Ka-band')
#ax2.set_title('W-band')
ax1.xaxis.set_major_formatter(xfmt)
ax2.xaxis.set_major_formatter(xfmt)
ax1.grid()
ax2.grid()
f.tight_layout(pad=0)
f.savefig(plotpath+'SW'+'.png', dpi=dpi, bbox_inches='tight')


plt.close('all')
