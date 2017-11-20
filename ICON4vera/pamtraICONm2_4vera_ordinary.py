from __future__ import division
from sys import argv, path
path.append('/home/dori/lib/python/')
import os
os.environ['PAMTRA_DATADIR'] = '/net/sever/mech/pamtra/data/'

import pyPamtra
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.dates as md
import time
from datetime import datetime

output_Z = 'output.png'
output_nc = 'output.nc'

#########################################################################
# PATHS
#########################################################################

# Directory where all the .nc meteograms are stored
ICON_folder = '/data/inscape/icon/experiments/nyalesund/newicon-2017-06-23-albedo/'

# Folder where your descriptor files are stored (you can use mine, Mario's or the default pamtra)
#descriptor_folder = '/home/mech/workspace/pamtra/descriptorfiles/'
descriptor_folder = '/home/dori/pamtra/descriptorfiles/'

#########################################################################
# FILES
#########################################################################

# Meteogram
ICON_filename = ICON_folder+'METEOGRAM_patch004_awipev.nc'

script, ICON_filename, output_nc, output_Z = argv


# Descriptor file for hydrometeors (Scattering models, m(D), v(D))
#descriptor_filename = 'descriptor_file_2m_liudb.txt'
descriptor_filename = 'descriptor_file_2m_ssrg.txt'

#########################################################################
# INIT
#########################################################################
pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptor_folder+descriptor_filename)

# SETTINGS
pam.nmlSet['active'] = True
pam.nmlSet['passive'] = False # Passive is time consuming
pam.set['verbose'] = 0 # set verbosity levels
pam.set['pyVerbose'] = 1 # change to 0 if you do not want to see job progress number

#########################################################################
# READ DATA
ICON_file = Dataset(ICON_filename, mode='r')
#########################################################################

vals = ICON_file.variables
Nh = len(vals['height_2'])
Nt = len(vals['time'])
H = np.tile(vals['height_2'][:],(Nt,1))
tt = np.tile(vals['time'][:],(Nh,1)).T
rain = vals['RAIN_GSP'][:] + vals['SNOW_GSP'][:]
print H.shape, tt.shape, rain.shape

# PLOT microphysic quantities (if you want)
#val = vals['QNR'][:]#vals['QNI'][:]+vals['QNS'][:]+vals['QNR'][:]+vals['QNG'][:]+vals['QNH'][:]+vals['QNC'][:] 
#plt.pcolormesh(tt,H,val)
#plt.ylim([0,3000])
#plt.colorbar()
#plt.show()

#timeidx = np.arange(500,600)
timeidx = np.arange(0,Nt)
pamData = dict() # empty dictionary to store pamtra Data

times   = vals['time'] # seconds since SOMETHING proplectic gregorian
pamData['hgt'] = np.tile(np.flip(vals['height_2'],0),(len(timeidx),1)) # heights at which fields are defined
units=times.units.split('since')[0]
basetime=pd.to_datetime(times.units.split('since')[-1])
dtimes = pd.to_timedelta(times[timeidx],unit=str(units)[0]) # TODO better than this...
pamData['timestamp'] = (basetime + dtimes).astype(np.int64)

print H.shape, tt.shape, rain.shape, pamData['timestamp'].shape
pamData['press']    = np.flip(vals['P'][timeidx],1)    # pressure
pamData['temp']     = np.flip(vals['T'][timeidx],1)    # temperature
pamData['relhum'] = np.flip(vals['REL_HUM'][timeidx],1)

# Read hydrometeors content
hydro_cmpl = np.zeros((len(timeidx),150,6))

hydro_cmpl[:,:,0] = np.flip(vals['QC'][timeidx],1)   # specific cloud water content
hydro_cmpl[:,:,1] = np.flip(vals['QI'][timeidx],1)   # specific cloud ice content
hydro_cmpl[:,:,2] = np.flip(vals['QR'][timeidx],1)   # rain mixing ratio
hydro_cmpl[:,:,3] = np.flip(vals['QS'][timeidx],1)   # snow mixing ratio
hydro_cmpl[:,:,4] = np.flip(vals['QG'][timeidx],1)   # graupel mixing ratio
hydro_cmpl[:,:,5] = np.flip(vals['QH'][timeidx],1)   # graupel mixing ratio

# Read hydrometeors number concentration
hydro_num_cmpl = np.zeros((len(timeidx),150,6))
hydro_num_cmpl[:,:,0] = np.flip(vals['QNC'][timeidx],1)  # number concentration of cloud water
hydro_num_cmpl[:,:,1] = np.flip(vals['QNI'][timeidx],1)  # number concentration ice
hydro_num_cmpl[:,:,2] = np.flip(vals['QNR'][timeidx],1)  # number concentration droplets
hydro_num_cmpl[:,:,3] = np.flip(vals['QNS'][timeidx],1)  # number concentration snow
hydro_num_cmpl[:,:,4] = np.flip(vals['QNG'][timeidx],1)  # number concentration graupel
hydro_num_cmpl[:,:,5] = np.flip(vals['QNH'][timeidx],1)  # number concentration hail


pamData["hydro_q"] = hydro_cmpl
pamData["hydro_n"] = hydro_num_cmpl

pam.createProfile(**pamData)

#########################################################################
# RUN
#########################################################################
frequencies = [94]
cores = 4 # number of parallel cores
pam.runParallelPamtra(np.array(frequencies), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF(output_nc) # SAVE OUTPUT

#########################################################################
# PLOTTING :)
#########################################################################
ZeW = pam.r['Ze'][:,0,:,0,0,0]
ZeW = np.ma.masked_values(ZeW,-9999)#[:-1,:]
H = pam.p['hgt'][:,0,:]#[:-1,:]
tt=np.tile((basetime + dtimes),(150,1)).T#[:-1,:]

levels = np.arange(-25,35,1)
print(tt.shape,H.shape,ZeW.shape)
plt.subplot(1, 1, 1)
plt.pcolormesh(tt,H,ZeW,vmin=-25,vmax=35,cmap='jet')
plt.text(0.1,0.1,'W-band', transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=1),fontsize='x-large')
plt.colorbar(label="W band Ze [dBZ]")
plt.grid()
plt.ylim([0,10000])
plt.ylabel('height [m]')
plt.tight_layout()
plt.savefig(output_Z)
