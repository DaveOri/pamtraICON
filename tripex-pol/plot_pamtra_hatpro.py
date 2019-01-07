# -*- coding: utf-8 -*-
"""
Import from pamtra simulated Tb and plot in the same format as hatpro quicklooks

"""
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import argparse

parser =  argparse.ArgumentParser(description='do plots for QuickLookBrowser')
parser.add_argument('-d','--date', nargs=1,
                    help='gimme datestring in the format YYYYMMDD')
parser.print_help()
args = parser.parse_args()
datestr = args.date[0]

runFld = '/data/optimice/pamtra_runs/tripex-pol/data/'
plotFld = '/data/optimice/pamtra_runs/tripex-pol/plots/'

xfmt = md.DateFormatter('%H')
plt.close('all')

pamtrafile = netCDF4.Dataset(runFld+datestr+'hatpro.nc')
datavars = pamtrafile.variables
datadims = pamtrafile.dimensions
datetime = netCDF4.num2date(datavars['datatime'][:],
                            datavars['datatime'].units)
timestamp = netCDF4.date2num(datetime, 'seconds since 1970-01-01 00:00:00')
tb = datavars['tb'][:,0,1,31,:,0] # downwelling at 0 meters

def plot_one_frequency(ax, time, tb, frequency, noxtick=True):
    ax.plot(time, tb,'k')
    ax.plot(time, tb.mean()*np.ones(time.shape), 'k', ls=':')
    ax.text(0.9, 0.8, str(tb.mean())+' K' ,
                  horizontalalignment='center',
                  verticalalignment='center',
                  color='red',
                  transform=ax.transAxes)
    ax.text(0.1, 0.8, str(frequency)+' GHz' ,
                  horizontalalignment='center',
                  verticalalignment='center',
                  color='red',
                  transform=ax.transAxes)
    ax.set_xlim([min(time),max(time)])
    ax.xaxis.set_major_formatter(xfmt)
    if noxtick:
        ax.set_xticklabels([])

f, axs = plt.subplots(8, 2, sharex=False, figsize=(11,9),
                      gridspec_kw = {'hspace':0.1})
axs[0,0].plot(datetime, 0.0*timestamp)
axs[0,1].plot(datetime, 0.0*timestamp)

plot_one_frequency(axs[1,0], datetime, tb[:,0], datavars['frequency'][0])
plot_one_frequency(axs[2,0], datetime, tb[:,1], datavars['frequency'][1])
plot_one_frequency(axs[3,0], datetime, tb[:,2], datavars['frequency'][2])
plot_one_frequency(axs[4,0], datetime, tb[:,3], datavars['frequency'][3])
plot_one_frequency(axs[5,0], datetime, tb[:,4], datavars['frequency'][4])
plot_one_frequency(axs[6,0], datetime, tb[:,5], datavars['frequency'][5])
plot_one_frequency(axs[7,0], datetime, tb[:,6], datavars['frequency'][6],False)
plot_one_frequency(axs[1,1], datetime, tb[:,7], datavars['frequency'][7])
plot_one_frequency(axs[2,1], datetime, tb[:,8], datavars['frequency'][8])
plot_one_frequency(axs[3,1], datetime, tb[:,9], datavars['frequency'][9])
plot_one_frequency(axs[4,1], datetime, tb[:,10], datavars['frequency'][10])
plot_one_frequency(axs[5,1], datetime, tb[:,11], datavars['frequency'][11])
plot_one_frequency(axs[6,1], datetime, tb[:,12], datavars['frequency'][12])
plot_one_frequency(axs[7,1], datetime, tb[:,13], datavars['frequency'][13],False)
f.tight_layout()
f.savefig(plotFld+datestr+'hatpro'+'.png', dpi=200, bbox_inches='tight')
plt.close('all')

