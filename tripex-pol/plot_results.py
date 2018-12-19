# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 17:34:01 2017

@author: dori
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
plt.close('all')
import matplotlib.dates as md
import numpy as np
import pandas as pd
from glob import glob
import os

import argparse
parser =  argparse.ArgumentParser(description='do plots for QuickLookBrowser')
parser.add_argument('--date', nargs=1, help='gimme datestring in the format YYYYMMDD')
parser.print_help()

args = parser.parse_args()
print(args)

plt.close('all')
MDV = False

runFld = '/data/optimice/pamtra_runs/tripex-pol/data/'
plotFld = '/data/optimice/pamtra_runs/tripex-pol/plots/'

runs = ['all_hydro_mom.nc', 'no_snow_mom.nc', 'only_ice_mom.nc', 'only_liquid_mom.nc', 'only_snow_mom.nc', 'only_graupel_hail_mom.nc']
titles = ['all Hydrometeors', 'No Snow', 'Only Ice', 'Only liquid (cloud drops and rain)', 'only Snow', 'only Graupel and Hail']
runTitles=dict(zip(runs,titles))

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

versus = -1 # Top Down
versus =  1 # Bottom Up

xfmt = md.DateFormatter('%m-%d %H')
ylim=(0,12000)
xDataLim = -1
figsize31=(18,18)
figsize21=(18,12)

def run41day(datestr):
    for run in runs[:]:
        print(datestr+' '+run)
        # Open the netcdf results file
        runFile = runFld + datestr + run
        runDataset = Dataset(runFile)
        runVars = runDataset.variables

        # Extract Heights
        H = (runVars['height'][:,0,:])[:xDataLim,:] # Kilometers
        # Reshape times
        ttt = pd.to_datetime(runVars['datatime'][:,0],unit='s')
        tt = (np.tile(ttt,(H.shape[1],1)).T)[:xDataLim,:]
        print(tt.shape, H.shape)
        # Extract attenuation
        ax = 2.0*(runVars['Attenuation_Hydrometeors'][:,0,:,0,0] + 
            runVars['Attenuation_Atmosphere'][:,0,:,0])
        Ax = ax[:,::versus].cumsum(axis=1)[:,::versus][:xDataLim,:]
        au = 2.0*(runVars['Attenuation_Hydrometeors'][:,0,:,1,0] + 
            runVars['Attenuation_Atmosphere'][:,0,:,1])
        Au = au[:,::versus].cumsum(axis=1)[:,::versus][:xDataLim,:]
        aa = 2.0*(runVars['Attenuation_Hydrometeors'][:,0,:,2,0] + 
            runVars['Attenuation_Atmosphere'][:,0,:,2])
        Aa = aa[:,::versus].cumsum(axis=1)[:,::versus][:xDataLim,:]
        aw = 2.0*(runVars['Attenuation_Hydrometeors'][:,0,:,3,0] + 
            runVars['Attenuation_Atmosphere'][:,0,:,3])
        Aw = aw[:,::versus].cumsum(axis=1)[:,::versus][:xDataLim,:]
        ag = 2.0*(runVars['Attenuation_Hydrometeors'][:,0,:,4,0] + 
            runVars['Attenuation_Atmosphere'][:,0,:,4])
        Ag = ag[:,::versus].cumsum(axis=1)[:,::versus][:xDataLim,:]

        # Plot Attenuation
        f,((ax1,ax2,ax3)) = plt.subplots(3, 1, sharex=False, figsize=figsize31)
        plot_variable(tt,H,Ax,ax1,None,'height [km]','dB','X-band 2-way Attenuation',0,1,ylim=ylim)
        plot_variable(tt,H,Aa,ax2,None,'height [km]','dB','Ka-band 2-way Attenuation',0,5,ylim=ylim)
        plot_variable(tt,H,Aw,ax3,'time','height [km]','dB', 'W-band 2-way Attenuation',0,15,ylim=ylim)
        #f.suptitle(runTitles[run], weight='black',bbox=dict(facecolor='white'))
        ax1.set_title('X-band')
        ax2.set_title('Ka-band')
        ax3.set_title('W-band')
        ax1.xaxis.set_major_formatter(xfmt)
        ax2.xaxis.set_major_formatter(xfmt)
        ax3.xaxis.set_major_formatter(xfmt)
        ax1.grid(color='k')
        ax2.grid(color='k')
        ax3.grid(color='k')
        f.tight_layout(pad=0)
        f.savefig(plotFld+datestr+run[:-3]+'_attenuation'+'.png', dpi=200, bbox_inches='tight')
        
        # Extract Ze
        Zex = runVars['Ze'][:,0,:,0,0,0][:xDataLim,:]
        Zeu = runVars['Ze'][:,0,:,1,0,0][:xDataLim,:]
        Zea = runVars['Ze'][:,0,:,2,0,0][:xDataLim,:]
        Zew = runVars['Ze'][:,0,:,3,0,0][:xDataLim,:]
        Zeg = runVars['Ze'][:,0,:,4,0,0][:xDataLim,:]

        # Extract Ze
        MDVx = -runVars['Radar_MeanDopplerVel'][:,0,:,0,0,0][:xDataLim,:]
        MDVu = -runVars['Radar_MeanDopplerVel'][:,0,:,1,0,0][:xDataLim,:]
        MDVa = -runVars['Radar_MeanDopplerVel'][:,0,:,2,0,0][:xDataLim,:]
        MDVw = -runVars['Radar_MeanDopplerVel'][:,0,:,3,0,0][:xDataLim,:]
        MDVg = -runVars['Radar_MeanDopplerVel'][:,0,:,4,0,0][:xDataLim,:]

        # Extract Ze
        SWx = runVars['Radar_SpectrumWidth'][:,0,:,0,0,0][:xDataLim,:]
        SWu = runVars['Radar_SpectrumWidth'][:,0,:,1,0,0][:xDataLim,:]
        SWa = runVars['Radar_SpectrumWidth'][:,0,:,2,0,0][:xDataLim,:]
        SWw = runVars['Radar_SpectrumWidth'][:,0,:,3,0,0][:xDataLim,:]
        SWg = runVars['Radar_SpectrumWidth'][:,0,:,4,0,0][:xDataLim,:]

        # Plot Ze
        f,((ax1,ax2,ax3)) = plt.subplots(3,1,sharex=False,figsize=figsize31)
        plot_variable(tt,H,Zex,ax1,None,'height [km]','dBZ','X-band Ze',-35,25,ylim=ylim)
        plot_variable(tt,H,Zea,ax2,None,'height [km]','dBZ', 'Ka-band Ze',-35,25,ylim=ylim)
        plot_variable(tt,H,Zew,ax3,'time','height [km]','dBZ', 'W-band Ze',-35,25,ylim=ylim)
        #f.suptitle(runTitles[run], weight='black',bbox=dict(facecolor='white'))
        ax1.xaxis.set_major_formatter(xfmt)
        ax2.xaxis.set_major_formatter(xfmt)
        ax3.xaxis.set_major_formatter(xfmt)
        ax1.set_title('X-band')
        ax2.set_title('Ka-band')
        ax3.set_title('W-band')
        ax1.grid(color='k')
        ax2.grid(color='k')
        ax3.grid(color='k')
        f.tight_layout(pad=0)
        f.savefig(plotFld+datestr+run[:-3]+'_Ze'+'.png', dpi=200, bbox_inches='tight')

        # make DWRs and plot
        DWRxa = Zex-Zea
        DWRaw = Zea-Zew
        DWRag = Zea-Zeg
        DWRwg = Zew-Zeg
        f,((ax1,ax2)) = plt.subplots(2,1,sharex=False,figsize=figsize21)
        plot_variable(tt,H,DWRxa,ax1,None,'height [km]','dB','DWR$_{X Ka}$',-5,20, ylim=ylim,cmap='nipy_spectral')
        plot_variable(tt,H,DWRaw,ax2,'time','height [km]','dB','DWR$_{Ka W}$',-5,20, ylim=ylim,cmap='nipy_spectral')
        f.suptitle(runTitles[run], weight='black',bbox=dict(facecolor='white'))
        ax1.xaxis.set_major_formatter(xfmt)
        ax2.xaxis.set_major_formatter(xfmt)
        ax1.set_title('X-Ka')
        ax2.set_title('Ka-W')
        ax1.grid(color='k')
        ax2.grid(color='k')
        f.tight_layout(pad=0)
        f.savefig(plotFld+datestr+run[:-3]+'_DWRe'+'.png', dpi=200, bbox_inches='tight')

        # make attenuated Z and DWRs and respective plots
        Zx = Zex-Ax
        Zu = Zeu-Au
        Za = Zea-Aa
        Zw = Zew-Aw
        Zg = Zeg-Ag
        f,((ax1,ax2,ax3)) = plt.subplots(3,1,sharex=False,figsize=figsize31)
        plot_variable(tt,H,Zx,ax1,None,'height [km]','dBZ','X-band Z attenuated',-35,25,ylim=ylim)
        plot_variable(tt,H,Za,ax2,None,'height [km]','dBZ','Ka-band Z attenuated',-35,25,ylim=ylim)
        plot_variable(tt,H,Zw,ax3,'time','height [km]','dBZ', 'W-band Z attenuated',-35,25,ylim=ylim)
        #f.suptitle(runTitles[run], weight='black',bbox=dict(facecolor='white'))
        ax1.set_title('X-band')
        ax2.set_title('Ka-band')
        ax3.set_title('W-band')
        ax1.xaxis.set_major_formatter(xfmt)
        ax2.xaxis.set_major_formatter(xfmt)
        ax3.xaxis.set_major_formatter(xfmt)
        ax1.grid(color='k')
        ax2.grid(color='k')
        ax3.grid(color='k')
        f.tight_layout(pad=0)
        f.savefig(plotFld+datestr+run[:-3]+'_Zattenuated'+'.png', dpi=200, bbox_inches='tight')

        DWRxa = Zx-Za
        DWRaw = Za-Zw
        DWRag = Za-Zg
        DWRwg = Zw-Zg
        f,((ax1,ax2)) = plt.subplots(2,1,sharex=False,figsize=figsize21)
        plot_variable(tt,H,DWRxa,ax1,None,'height [km]','dB','DWR$_{X Ka}$ attenuated',-5,20,ylim=ylim,cmap='nipy_spectral')    
        plot_variable(tt,H,DWRaw,ax2,'time','height [km]','dB','DWR$_{Ka W}$ attenuated',-5,20,ylim=ylim,cmap='nipy_spectral')
        #f.suptitle(runTitles[run], weight='black',bbox=dict(facecolor='white'))
        ax1.set_title('X-Ka')
        ax2.set_title('Ka-W')
        ax1.xaxis.set_major_formatter(xfmt)
        ax2.xaxis.set_major_formatter(xfmt)
        ax1.grid(color='k')
        ax2.grid(color='k')
        f.tight_layout(pad=0)
        f.savefig(plotFld+datestr+run[:-3]+'_DWRattenuated'+'.png', dpi=200, bbox_inches='tight')
        
        # Plot mean doppler velocity
        f,((ax1,ax2,ax3)) = plt.subplots(3,1,sharex=False,figsize=figsize31)
        plot_variable(tt,H,MDVx,ax1,None,  'height [km]','m/s','Ku-band MDV',-3,0,ylim=ylim)
        plot_variable(tt,H,MDVa,ax2,None,  'height [km]','m/s','Ka-band MDV',-3,0,ylim=ylim)
        plot_variable(tt,H,MDVw,ax3,'time','height [km]','m/s', 'W-band MDV',-3,0,ylim=ylim)
        #f.suptitle(runTitles[run], weight='black',bbox=dict(facecolor='white'))
        ax1.set_title('X-band')
        ax2.set_title('Ka-band')
        ax3.set_title('W-band')
        ax1.xaxis.set_major_formatter(xfmt)
        ax2.xaxis.set_major_formatter(xfmt)
        ax3.xaxis.set_major_formatter(xfmt)
        ax1.grid(color='k')
        ax2.grid(color='k')
        ax3.grid(color='k')
        f.tight_layout(pad=0)
        f.savefig(plotFld+datestr+run[:-3]+'_MDV'+'.png', dpi=200, bbox_inches='tight')
      
        f,((ax1,ax2,ax3)) = plt.subplots(3,1,sharex=False,figsize=figsize31)
        plot_variable(tt,H,SWx,ax1,None,  'height [km]','m/s','Ku-band SW',0,1,ylim=ylim)
        plot_variable(tt,H,SWa,ax2,None,  'height [km]','m/s','Ka-band SW',0,1,ylim=ylim)
        plot_variable(tt,H,SWw,ax3,'time','height [km]','m/s', 'W-band SW',0,1,ylim=ylim)
        #f.suptitle(runTitles[run], weight='black',bbox=dict(facecolor='white'))
        ax1.set_title('X-band')
        ax2.set_title('Ka-band')
        ax3.set_title('W-band')
        ax1.xaxis.set_major_formatter(xfmt)
        ax2.xaxis.set_major_formatter(xfmt)
        ax3.xaxis.set_major_formatter(xfmt)
        ax1.grid(color='k')
        ax2.grid(color='k')
        ax3.grid(color='k')
        f.tight_layout(pad=0)
        f.savefig(plotFld+datestr+run[:-3]+'_SW'+'.png', dpi=200, bbox_inches='tight')

        # Plot dual doppler velocity
        DDWxa = MDVx-MDVa
        DDWaw = MDVa-MDVw
        DDWag = MDVa-MDVg
        DDWwg = MDVw-MDVg
        f,((ax1,ax2)) = plt.subplots(2,1,sharex=False,figsize=figsize21)
        plot_variable(tt,H,DDWxa,ax1,None,'height [km]','m/s','DDV$_{X Ka}$',-0.3,0.3,ylim=ylim,cmap='nipy_spectral')
        plot_variable(tt,H,DDWaw,ax2,'time','height [km]','m/s','DDV$_{Ka W}$',-0.3,0.3,ylim=ylim,cmap='nipy_spectral')
        #f.suptitle(runTitles[run], weight='black',bbox=dict(facecolor='white'))
        ax1.set_title('X-Ka')
        ax2.set_title('Ka-W')
        ax1.xaxis.set_major_formatter(xfmt)
        ax2.xaxis.set_major_formatter(xfmt)
        ax1.grid(color='k')
        ax2.grid(color='k')
        f.tight_layout(pad=0)
        f.savefig(plotFld+datestr+run[:-3]+'_DDV'+'.png', dpi=200, bbox_inches='tight')

        # Plot dual spectral width
        DSWxa = SWx-SWa
        DSWaw = SWa-SWw
        DSWag = SWa-SWg
        DSWwg = SWw-SWg
        f,((ax1,ax2)) = plt.subplots(2,1,sharex=False,figsize=figsize21)
        plot_variable(tt,H,DSWxa,ax1,None,'height [km]','m/s','DSW$_{X Ka}$',-0.3,0.3,ylim=ylim,cmap='nipy_spectral')
        plot_variable(tt,H,DSWaw,ax2,'time','height [km]','m/s','DSW$_{Ka W}$',-0.3,0.3,ylim=ylim,cmap='nipy_spectral')
        #f.suptitle(runTitles[run], weight='black',bbox=dict(facecolor='white'))
        ax1.set_title('X-Ka')
        ax2.set_title('Ka-W')
        ax1.grid(color='k')
        ax2.grid(color='k')
        ax1.xaxis.set_major_formatter(xfmt)
        ax2.xaxis.set_major_formatter(xfmt)
        f.tight_layout(pad=0)
        f.savefig(plotFld+datestr+run[:-3]+'_DSW'+'.png', dpi=200, bbox_inches='tight')

        plt.close('all')

if (args.date[0] == 'all'):
    alldates = [os.path.basename(x)[:8] for x in glob('/data/optimice/pamtra_runs/tripex-pol/data/*all_hydro_mom.nc')]
    print('Running over all dates',alldates)
    for dd in alldates:
        run41day(dd)
else:
    run41day(datestr = args.date[0])
