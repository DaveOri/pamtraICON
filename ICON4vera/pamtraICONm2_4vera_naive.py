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
ICON_folder = '/data/inscape/icon/experiments/nyalesund/newicon-2017-06-23-albedo/postproc/'

# Folder where your descriptor files are stored (you can use mine, Mario's or the default pamtra)
#descriptor_folder = '/home/mech/workspace/pamtra/descriptorfiles/'
descriptor_folder = '/home/dori/pamtra/descriptorfiles/'
descriptor_folder = '/home/dori/descriptorfiles/'

#########################################################################
# FILES
#########################################################################

# Meteogram
ICON_filename = ICON_folder + '599_leg1nya.nc'

script, ICON_filename, output_nc, output_Z = argv

# Descriptor file for hydrometeors (Scattering models, m(D), v(D))
#descriptor_filename = 'descriptor_file_2m_liudb.txt'
descriptor_filename = 'descriptor_file_2m_ssrg.txt'
descriptor_filename = 'descriptor_file_2m_ssrgNEW.txt'

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
Nh = len(vals['height'])
Nt = len(vals['clat'])
H = vals['height_2d'][:].T
tt = np.tile(vals['clat'][:],(Nh,1)).T
print H.shape, tt.shape

# PLOT microphysic quantities (if you want)
#val = vals['QNR'][:]#vals['QNI'][:]+vals['QNS'][:]+vals['QNR'][:]+vals['QNG'][:]+vals['QNH'][:]+vals['QNC'][:] 
#plt.pcolormesh(tt,H,val)
#plt.ylim([0,3000])
#plt.colorbar()
#plt.show()

#timeidx = np.arange(500,600)
timeidx = np.arange(0,Nt)
pamData = dict() # empty dictionary to store pamtra Data

clat   = vals['clat']
clon   = vals['clon']
pamData['hgt'] = np.flip(H,1) #np.tile(np.flip(vals['height'],0),(len(timeidx),1)) # heights at which fields are defined

pamData['press']  = np.flip(vals['P'][:].T,1)    # pressure
pamData['temp']   = np.flip(vals['T'][:].T,1)    # temperature
pamData['relhum'] = np.flip(vals['REL_HUM'][:].T,1)

# Read hydrometeors content
hydro_cmpl = np.zeros((Nt,Nh,6))
#hydro_cmpl[:,:,#] = vals['QV'][:]   # specific humidity
hydro_cmpl[:,:,0] = np.flip(vals['QC'][:].T,1) # specific cloud water content
hydro_cmpl[:,:,1] = np.flip(vals['QI'][:].T,1) # specific cloud ice content
hydro_cmpl[:,:,2] = np.flip(vals['QR'][:].filled(fill_value=0.0).T,1) # rain mixing ratio
hydro_cmpl[:,:,3] = np.flip(vals['QS'][:].T,1) # snow mixing ratio
hydro_cmpl[:,:,4] = np.flip(vals['QG'][:].T,1) # graupel mixing ratio
hydro_cmpl[:,:,5] = np.flip(vals['QH'][:].T,1)   # graupel mixing ratio # TODO report probably error encoding long name ...  should be hail mixing ratio

# Read hydrometeors number concentration
hydro_num_cmpl = np.zeros((len(timeidx),Nh,6))
hydro_num_cmpl[:,:,0] = np.flip(vals['QNC'][:].T,1)  # number concentration of cloud water
hydro_num_cmpl[:,:,1] = np.flip(vals['QNI'][:].T,1)  # number concentration ice
hydro_num_cmpl[:,:,2] = np.flip(vals['QNR'][:].filled(fill_value=0.0).T,1)  # number concentration droplets
hydro_num_cmpl[:,:,3] = np.flip(vals['QNS'][:].T,1)  # number concentration snow
hydro_num_cmpl[:,:,4] = np.flip(vals['QNG'][:].T,1)  # number concentration graupel
hydro_num_cmpl[:,:,5] = np.flip(vals['QNH'][:].T,1)  # number concentration hail

pamData["hydro_q"] = hydro_cmpl
pamData["hydro_n"] = hydro_num_cmpl

pam.createProfile(**pamData)

#########################################################################
# RUN
#########################################################################
frequencies = [9.6,13.6,35.6,94,220]
cores = 1 # number of parallel cores
pam.runParallelPamtra(np.array(frequencies), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF(output_nc) # SAVE OUTPUT

#########################################################################
# PLOTTING :)
#########################################################################
ZeW = pam.r['Ze'][:,0,:,3,0,0]
#ZeW = pam.r['Ze'][:,0,:,2,0,0]
ZeW = np.ma.masked_values(ZeW,-9999)

H = pam.p['hgt'][:,0,:]
tt= clat
tt=np.tile(tt[:].reshape((len(tt),1)),(1,H.shape[1]))

levels = np.arange(-25,35,1)
plt.subplot(1, 1, 1)
#plt.contourf(tt,H,ZeX,levels,cmap='jet')
plt.pcolormesh(tt,H,ZeW,vmin=-25,vmax=35,cmap='jet')
plt.text(0.1,0.9,'W-band', transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=1),fontsize='x-large')
plt.colorbar(label="W band Ze [dBZ]")
plt.grid()
plt.ylim([0,10000])
plt.ylabel('height [m]')
plt.title(descriptor_filename.split('_')[-1])
plt.savefig(output_Z)
