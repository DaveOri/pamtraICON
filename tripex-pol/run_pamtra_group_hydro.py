from __future__ import division
import pyPamtra
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.dates as md
import time
from datetime import datetime
from sys import argv, path
import gc

import argparse
parser =  argparse.ArgumentParser(description='do plots for QuickLookBrowser')
parser.add_argument('--date', nargs=1, help='gimme datestring in the format YYYYMMDD')
parser.print_help()
args = parser.parse_args()
print(args)

cores = 4 # number of parallel cores

#########################################################################
# PATHS
#########################################################################

# Directory where all the .nc meteograms are stored
datestr = args.date[0]
ICON_folder = '/data/inscape/icon/experiments/juelich/testbed/testbed_' + datestr + '/'

# Folder where your descriptor files are stored (you can use mine, Mario's or the default pamtra)
descriptor_folder = '/home/dori/descriptorfiles/'

#########################################################################
# FILES
#########################################################################

# Meteogram
ICON_filename = ICON_folder + 'METEOGRAM_patch001_' + datestr + '_joyce.nc'

# Descriptor file for hydrometeors (Scattering models, m(D), v(D))
#descriptor_filename = 'descriptor_file_2m_ssrgNEW.txt'
descriptor_filename = 'descriptor_file_2m_ssrgNEWpowerLaw.txt'

# Directory of pluvios for precipitation comparison
#plufile = '/data/data_hatpro/jue/data/pluvio/netcdf/1511/pluvio2_jue_20151124.nc'
#pluvio = Dataset(plufile, 'r')
#pluvars = pluvio.variables
#var = 'total_accum_NRT'
##var = 'fill_level_NRT'
#TotAccNRT = pluvars[var]
#PluvioTime = pluvars['time']
#units=PluvioTime.units.split('since')[0]
#btPluvio = pd.to_datetime(PluvioTime.units.split('since')[-1])
#dtPluvio = pd.to_timedelta([i if i > 0 else np.nan for i in PluvioTime[:]],unit=str(units)[0]) # TODO better than this...
#PluvioDateTime = (btPluvio + dtPluvio)#.astype(np.int64)
#TotAcc= TotAccNRT[:] - TotAccNRT[0]


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
pam.nmlSet["radar_mode"] = "moments"
pam.nmlSet["radar_pnoise0"] = -60.0

#########################################################################
# READ DATA
print ICON_filename
ICON_file = Dataset(ICON_filename, mode='r')
#########################################################################
# Akio ICON files give all of the varaibles in a multidimensional array
# ICON.variables["values"] (5-D array)
# ICON files from Tripex are stored in a much more structured way
# LET'S DO IT !!!

vals = ICON_file.variables
Nh = len(vals['height_2'])
Nt = len(vals['time'])
H = np.tile(vals['height_2'][:],(Nt,1))
tt = np.tile(vals['time'][:],(Nh,1)).T
rain = vals['RAIN_GSP'][:] + vals['SNOW_GSP'][:] # rapid check against pluvio of rain accumulation at the ground
print H.shape, tt.shape, rain.shape

#timeidx = np.arange(4500,4504) 
timeidx = np.arange(0,Nt)
pamData = dict() # empty dictionary to store pamtra Data

times   = vals['time'] # seconds since 2015-11-24 02:00:03 proplectic gregorian
pamData['hgt'] = np.tile(np.flip(vals['height_2'],0),(len(timeidx),1)) # heights at which fields are defined
units=times.units.split('since')[0]
basetime=pd.to_datetime(times.units.split('since')[-1])
dtimes = pd.to_timedelta(times[timeidx],unit=str(units)[0]) # TODO better than this...
pamData['timestamp'] = (basetime + dtimes).astype(np.int64)*1e-9
#pamData['unixtime'] = (basetime + dtimes).astype(np.int64)*1.0e-9

print H.shape, tt.shape, rain.shape, #pamData['timestamp'].shape
pamData['press']    = np.flip(vals['P'][timeidx],1)    # pressure 
pamData['temp']     = np.flip(vals['T'][timeidx],1)    # temperature
#pamData['wind_u']  = vals['U'][:]    # zonal wind speed
#pamData['wind_v']  = vals['V'][:]    # meridional wind speed
wind_w = np.flip(vals['W'][timeidx],1)    # vertical wind speed
pamData['wind_w'] = 0.5*(wind_w[:,:-1]+wind_w[:,1:])
pamData['relhum'] = np.flip(vals['REL_HUM'][timeidx],1)

# Read hydrometeors content
hydro_cmpl = np.zeros((len(timeidx),Nh,6))
#hydro_cmpl[:,:,#] = vals['QV'][:]   # specific humidity
hydro_cmpl[:,:,0] = 0.0#np.flip(vals['QC'][timeidx],1)   # specific cloud water content
hydro_cmpl[:,:,1] = 0.0#np.flip(vals['QI'][timeidx],1)   # specific cloud ice content
hydro_cmpl[:,:,2] = 0.0#np.flip(vals['QR'][timeidx],1)   # rain mixing ratio
hydro_cmpl[:,:,3] = np.flip(vals['QS'][timeidx],1)   # snow mixing ratio
hydro_cmpl[:,:,4] = 0.0#np.flip(vals['QG'][timeidx],1)   # graupel mixing ratio
hydro_cmpl[:,:,5] = 0.0#np.flip(vals['QH'][timeidx],1)   # graupel mixing ratio # TODO report probably error encoding long name ...  should be hail mixing ratio

# Read hydrometeors number concentration
hydro_num_cmpl = np.zeros((len(timeidx),Nh,6))
hydro_num_cmpl[:,:,0] = 0.0#np.flip(vals['QNC'][timeidx],1)  # number concentration of cloud water
hydro_num_cmpl[:,:,1] = 0.0#np.flip(vals['QNI'][timeidx],1)  # number concentration ice
hydro_num_cmpl[:,:,2] = 0.0#np.flip(vals['QNR'][timeidx],1)  # number concentration droplets
hydro_num_cmpl[:,:,3] = np.flip(vals['QNS'][timeidx],1)  # number concentration snow
hydro_num_cmpl[:,:,4] = 0.0#np.flip(vals['QNG'][timeidx],1)  # number concentration graupel
hydro_num_cmpl[:,:,5] = 0.0#np.flip(vals['QNH'][timeidx],1)  # number concentration hail 


pamData["hydro_q"] = hydro_cmpl
pamData["hydro_n"] = hydro_num_cmpl

#pamData["lfrac"] = np.array([1]) ####
#pamData["groundtemp"] = pamData['temp'][149]
#pamData['phalf'] = vals['PHALF'][:]# pressure on the half levels

pam.createProfile(**pamData)
#print(pam.p['unixtime'])
#print(basetime + dtimes)
#########################################################################
# RUN
#########################################################################
frequencies = [9.6,13.6,35.6,94,220]
pam.runParallelPamtra(np.array(frequencies), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF('/data/optimice/pamtra_runs/tripex-pol/data/'+datestr+'only_snow_mom.nc') # SAVE OUTPUT
gc.collect()

################################################################################
pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptor_folder+descriptor_filename)

# SETTINGS
pam.nmlSet['active'] = True
pam.nmlSet['passive'] = False # Passive is time consuming
pam.set['verbose'] = 0 # set verbosity levels
pam.set['pyVerbose'] = 1 # change to 0 if you do not want to see job progress number
pam.nmlSet["radar_mode"] = "moments"
pam.nmlSet["radar_pnoise0"] = -60.0

ICON_file = Dataset(ICON_filename, mode='r')

vals = ICON_file.variables
Nh = len(vals['height_2'])
Nt = len(vals['time'])
H = np.tile(vals['height_2'][:],(Nt,1))
tt = np.tile(vals['time'][:],(Nh,1)).T
rain = vals['RAIN_GSP'][:] + vals['SNOW_GSP'][:] # rapid check against pluvio of rain accumulation at the ground
print H.shape, tt.shape, rain.shape

# Read hydrometeors content
hydro_cmpl = np.zeros((len(timeidx),Nh,6))
hydro_cmpl[:,:,0] = np.flip(vals['QC'][timeidx],1)   # specific cloud water content
hydro_cmpl[:,:,1] = np.flip(vals['QI'][timeidx],1)   # specific cloud ice content
hydro_cmpl[:,:,2] = np.flip(vals['QR'][timeidx],1)   # rain mixing ratio
hydro_cmpl[:,:,3] = 0.0*np.flip(vals['QS'][timeidx],1)   # snow mixing ratio
hydro_cmpl[:,:,4] = np.flip(vals['QG'][timeidx],1)   # graupel mixing ratio
hydro_cmpl[:,:,5] = np.flip(vals['QH'][timeidx],1)   # graupel mixing ratio # TODO report probably error encoding long name ...  should be hail mixing ratio

# Read hydrometeors number concentration
hydro_num_cmpl = np.zeros((len(timeidx),Nh,6))
hydro_num_cmpl[:,:,0] = np.flip(vals['QNC'][timeidx],1)  # number concentration of cloud water
hydro_num_cmpl[:,:,1] = np.flip(vals['QNI'][timeidx],1)  # number concentration ice
hydro_num_cmpl[:,:,2] = np.flip(vals['QNR'][timeidx],1)  # number concentration droplets
hydro_num_cmpl[:,:,3] = 0.0*np.flip(vals['QNS'][timeidx],1)  # number concentration snow
hydro_num_cmpl[:,:,4] = np.flip(vals['QNG'][timeidx],1)  # number concentration graupel
hydro_num_cmpl[:,:,5] = np.flip(vals['QNH'][timeidx],1)  # number concentration hail 


pamData["hydro_q"] = hydro_cmpl
pamData["hydro_n"] = hydro_num_cmpl
pam.createProfile(**pamData)

frequencies = [9.6,13.6,35.6,94,220]
pam.runParallelPamtra(np.array(frequencies), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF('/data/optimice/pamtra_runs/tripex-pol/data/'+datestr+'no_snow_mom.nc') # SAVE OUTPUT
gc.collect()

################################################################################
pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptor_folder+descriptor_filename)

# SETTINGS
pam.nmlSet['active'] = True
pam.nmlSet['passive'] = False # Passive is time consuming
pam.set['verbose'] = 0 # set verbosity levels
pam.set['pyVerbose'] = 1 # change to 0 if you do not want to see job progress number
pam.nmlSet["radar_mode"] = "moments"
pam.nmlSet["radar_pnoise0"] = -60.0

ICON_file = Dataset(ICON_filename, mode='r')

vals = ICON_file.variables
Nh = len(vals['height_2'])
Nt = len(vals['time'])
H = np.tile(vals['height_2'][:],(Nt,1))
tt = np.tile(vals['time'][:],(Nh,1)).T
rain = vals['RAIN_GSP'][:] + vals['SNOW_GSP'][:] # rapid check against pluvio of rain accumulation at the ground
print H.shape, tt.shape, rain.shape

# Read hydrometeors content
hydro_cmpl = np.zeros((len(timeidx),Nh,6))
hydro_cmpl[:,:,0] = np.flip(vals['QC'][timeidx],1)   # specific cloud water content
hydro_cmpl[:,:,1] = 0.0*np.flip(vals['QI'][timeidx],1)   # specific cloud ice content
hydro_cmpl[:,:,2] = np.flip(vals['QR'][timeidx],1)   # rain mixing ratio
hydro_cmpl[:,:,3] = 0.0*np.flip(vals['QS'][timeidx],1)   # snow mixing ratio
hydro_cmpl[:,:,4] = 0.0*np.flip(vals['QG'][timeidx],1)   # graupel mixing ratio
hydro_cmpl[:,:,5] = 0.0*np.flip(vals['QH'][timeidx],1)   # graupel mixing ratio # TODO report probably error encoding long name ...  should be hail mixing ratio

# Read hydrometeors number concentration
hydro_num_cmpl = np.zeros((len(timeidx),Nh,6))
hydro_num_cmpl[:,:,0] = np.flip(vals['QNC'][timeidx],1)  # number concentration of cloud water
hydro_num_cmpl[:,:,1] = 0.0*np.flip(vals['QNI'][timeidx],1)  # number concentration ice
hydro_num_cmpl[:,:,2] = np.flip(vals['QNR'][timeidx],1)  # number concentration droplets
hydro_num_cmpl[:,:,3] = 0.0*np.flip(vals['QNS'][timeidx],1)  # number concentration snow
hydro_num_cmpl[:,:,4] = 0.0*np.flip(vals['QNG'][timeidx],1)  # number concentration graupel
hydro_num_cmpl[:,:,5] = 0.0*np.flip(vals['QNH'][timeidx],1)  # number concentration hail 


pamData["hydro_q"] = hydro_cmpl
pamData["hydro_n"] = hydro_num_cmpl
pam.createProfile(**pamData)

frequencies = [9.6,13.6,35.6,94,220]
pam.runParallelPamtra(np.array(frequencies), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF('/data/optimice/pamtra_runs/tripex-pol/data/'+datestr+'only_liquid_mom.nc') # SAVE OUTPUT
gc.collect()

################################################################################
pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptor_folder+descriptor_filename)

# SETTINGS
pam.nmlSet['active'] = True
pam.nmlSet['passive'] = False # Passive is time consuming
pam.set['verbose'] = 0 # set verbosity levels
pam.set['pyVerbose'] = 1 # change to 0 if you do not want to see job progress number
pam.nmlSet["radar_mode"] = "moments"
pam.nmlSet["radar_pnoise0"] = -60.0

ICON_file = Dataset(ICON_filename, mode='r')

vals = ICON_file.variables
Nh = len(vals['height_2'])
Nt = len(vals['time'])
H = np.tile(vals['height_2'][:],(Nt,1))
tt = np.tile(vals['time'][:],(Nh,1)).T
rain = vals['RAIN_GSP'][:] + vals['SNOW_GSP'][:] # rapid check against pluvio of rain accumulation at the ground
print H.shape, tt.shape, rain.shape

# Read hydrometeors content
hydro_cmpl = np.zeros((len(timeidx),Nh,6))
hydro_cmpl[:,:,0] = 0.0*np.flip(vals['QC'][timeidx],1)   # specific cloud water content
hydro_cmpl[:,:,1] = np.flip(vals['QI'][timeidx],1)   # specific cloud ice content
hydro_cmpl[:,:,2] = 0.0*np.flip(vals['QR'][timeidx],1)   # rain mixing ratio
hydro_cmpl[:,:,3] = 0.0*np.flip(vals['QS'][timeidx],1)   # snow mixing ratio
hydro_cmpl[:,:,4] = 0.0*np.flip(vals['QG'][timeidx],1)   # graupel mixing ratio
hydro_cmpl[:,:,5] = 0.0*np.flip(vals['QH'][timeidx],1)   # graupel mixing ratio # TODO report probably error encoding long name ...  should be hail mixing ratio

# Read hydrometeors number concentration
hydro_num_cmpl = np.zeros((len(timeidx),Nh,6))
hydro_num_cmpl[:,:,0] = 0.0*np.flip(vals['QNC'][timeidx],1)  # number concentration of cloud water
hydro_num_cmpl[:,:,1] = np.flip(vals['QNI'][timeidx],1)  # number concentration ice
hydro_num_cmpl[:,:,2] = 0.0*np.flip(vals['QNR'][timeidx],1)  # number concentration droplets
hydro_num_cmpl[:,:,3] = 0.0*np.flip(vals['QNS'][timeidx],1)  # number concentration snow
hydro_num_cmpl[:,:,4] = 0.0*np.flip(vals['QNG'][timeidx],1)  # number concentration graupel
hydro_num_cmpl[:,:,5] = 0.0*np.flip(vals['QNH'][timeidx],1)  # number concentration hail 


pamData["hydro_q"] = hydro_cmpl
pamData["hydro_n"] = hydro_num_cmpl
pam.createProfile(**pamData)

frequencies = [9.6,13.6,35.6,94,220]
pam.runParallelPamtra(np.array(frequencies), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF('/data/optimice/pamtra_runs/tripex-pol/data/'+datestr+'only_ice_mom.nc') # SAVE OUTPUT
gc.collect()

################################################################################
pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptor_folder+descriptor_filename)

# SETTINGS
pam.nmlSet['active'] = True
pam.nmlSet['passive'] = False # Passive is time consuming
pam.set['verbose'] = 0 # set verbosity levels
pam.set['pyVerbose'] = 1 # change to 0 if you do not want to see job progress number
pam.nmlSet["radar_mode"] = "moments"
pam.nmlSet["radar_pnoise0"] = -60.0

ICON_file = Dataset(ICON_filename, mode='r')

vals = ICON_file.variables
Nh = len(vals['height_2'])
Nt = len(vals['time'])
H = np.tile(vals['height_2'][:],(Nt,1))
tt = np.tile(vals['time'][:],(Nh,1)).T
rain = vals['RAIN_GSP'][:] + vals['SNOW_GSP'][:] # rapid check against pluvio of rain accumulation at the ground
print H.shape, tt.shape, rain.shape

# Read hydrometeors content
hydro_cmpl = np.zeros((len(timeidx),Nh,6))
hydro_cmpl[:,:,0] = 0.0*np.flip(vals['QC'][timeidx],1)   # specific cloud water content
hydro_cmpl[:,:,1] = 0.0*np.flip(vals['QI'][timeidx],1)   # specific cloud ice content
hydro_cmpl[:,:,2] = 0.0*np.flip(vals['QR'][timeidx],1)   # rain mixing ratio
hydro_cmpl[:,:,3] = 0.0*np.flip(vals['QS'][timeidx],1)   # snow mixing ratio
hydro_cmpl[:,:,4] = np.flip(vals['QG'][timeidx],1)   # graupel mixing ratio
hydro_cmpl[:,:,5] = np.flip(vals['QH'][timeidx],1)   # graupel mixing ratio # TODO report probably error encoding long name ...  should be hail mixing ratio

# Read hydrometeors number concentration
hydro_num_cmpl = np.zeros((len(timeidx),Nh,6))
hydro_num_cmpl[:,:,0] = 0.0*np.flip(vals['QNC'][timeidx],1)  # number concentration of cloud water
hydro_num_cmpl[:,:,1] = 0.0*np.flip(vals['QNI'][timeidx],1)  # number concentration ice
hydro_num_cmpl[:,:,2] = 0.0*np.flip(vals['QNR'][timeidx],1)  # number concentration droplets
hydro_num_cmpl[:,:,3] = 0.0*np.flip(vals['QNS'][timeidx],1)  # number concentration snow
hydro_num_cmpl[:,:,4] = np.flip(vals['QNG'][timeidx],1)  # number concentration graupel
hydro_num_cmpl[:,:,5] = np.flip(vals['QNH'][timeidx],1)  # number concentration hail 


pamData["hydro_q"] = hydro_cmpl
pamData["hydro_n"] = hydro_num_cmpl
pam.createProfile(**pamData)

frequencies = [9.6,13.6,35.6,94,220]
pam.runParallelPamtra(np.array(frequencies), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF('/data/optimice/pamtra_runs/tripex-pol/data/'+datestr+'only_graupel_hail_mom.nc') # SAVE OUTPUT
gc.collect()

#################################################################################

pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptor_folder+descriptor_filename)

# SETTINGS
pam.nmlSet['active'] = True
pam.nmlSet['passive'] = False # Passive is time consuming
pam.set['verbose'] = 0 # set verbosity levels
pam.set['pyVerbose'] = 1 # change to 0 if you do not want to see job progress number
pam.nmlSet["radar_mode"] = "moments"
pam.nmlSet["radar_pnoise0"] = -60.0

ICON_file = Dataset(ICON_filename, mode='r')

vals = ICON_file.variables
Nh = len(vals['height_2'])
Nt = len(vals['time'])
H = np.tile(vals['height_2'][:],(Nt,1))
tt = np.tile(vals['time'][:],(Nh,1)).T
rain = vals['RAIN_GSP'][:] + vals['SNOW_GSP'][:] # rapid check against pluvio of rain accumulation at the ground
print H.shape, tt.shape, rain.shape

# Read hydrometeors content
hydro_cmpl = np.zeros((len(timeidx),Nh,6))
hydro_cmpl[:,:,0] = np.flip(vals['QC'][timeidx],1)   # specific cloud water content
hydro_cmpl[:,:,1] = np.flip(vals['QI'][timeidx],1)   # specific cloud ice content
hydro_cmpl[:,:,2] = np.flip(vals['QR'][timeidx],1)   # rain mixing ratio
hydro_cmpl[:,:,3] = np.flip(vals['QS'][timeidx],1)   # snow mixing ratio
hydro_cmpl[:,:,4] = np.flip(vals['QG'][timeidx],1)   # graupel mixing ratio
hydro_cmpl[:,:,5] = np.flip(vals['QH'][timeidx],1)   # graupel mixing ratio # TODO report probably error encoding long name ...  should be hail mixing ratio

# Read hydrometeors number concentration
hydro_num_cmpl = np.zeros((len(timeidx),Nh,6))
hydro_num_cmpl[:,:,0] = np.flip(vals['QNC'][timeidx],1)  # number concentration of cloud water
hydro_num_cmpl[:,:,1] = np.flip(vals['QNI'][timeidx],1)  # number concentration ice
hydro_num_cmpl[:,:,2] = np.flip(vals['QNR'][timeidx],1)  # number concentration droplets
hydro_num_cmpl[:,:,3] = np.flip(vals['QNS'][timeidx],1)  # number concentration snow
hydro_num_cmpl[:,:,4] = np.flip(vals['QNG'][timeidx],1)  # number concentration graupel
hydro_num_cmpl[:,:,5] = np.flip(vals['QNH'][timeidx],1)  # number concentration hail 


pamData["hydro_q"] = hydro_cmpl
pamData["hydro_n"] = hydro_num_cmpl
pam.createProfile(**pamData)

frequencies = [9.6,13.6,35.6,94,220]
pam.runParallelPamtra(np.array(frequencies), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF('/data/optimice/pamtra_runs/tripex-pol/data/'+datestr+'all_hydro_mom.nc') # SAVE OUTPUT
gc.collect()
