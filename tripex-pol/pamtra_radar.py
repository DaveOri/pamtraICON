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

import argparse
parser =  argparse.ArgumentParser(description='do plots for QuickLookBrowser')
parser.add_argument('--date', nargs=1, help='gimme datestring in the format YYYYMMDD')
parser.print_help()
args = parser.parse_args()
print(args)

cores = 2 # number of parallel cores

#########################################################################
# PATHS
#########################################################################

# Directory where all the .nc meteograms are stored
datestr = args.date[0]
ICON_folder = '/data/inscape/icon/experiments/juelich/testbed/testbed_' + datestr + '/'

# Folder where your descriptor files are stored (you can use mine, Mario's or the default pamtra)
descriptor_folder = '/home/dori/descriptorfiles/'

#########################################################################
# INIT
#########################################################################
descriptor_filename = 'descriptor_file_2m_ssrgNEWpowerLaw.txt'
descriptorFile = np.array([
      #['hydro_name' 'as_ratio' 'liq_ice' 'rho_ms' 'a_ms' 'b_ms' 'alpha_as' 'beta_as' 'moment_in' 'nbin' 'dist_name' 'p_1' 'p_2' 'p_3' 'p_4' 'd_1' 'd_2' 'scat_name' 'vel_size_mod' 'canting']
       ('cwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,   2.0e-6,   8.0e-5,       'mie-sphere',     'khvorostyanov01_drops', -99.0),
       ('iwc_q', 1.0, -1, -99.0, 1.58783,  2.56,  0.684,   2.0, 13, 100, 'mgamma', -99.0, -99.0, 1.564, 0.8547, 1.744e-5, 9.369e-3, 'ss-rayleigh-gans', 'corPowerLaw_30.606_0.5533', -99.0),
       ('rwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,  0.00012,   8.2e-3,       'mie-sphere',     'khvorostyanov01_drops', -99.0),
       ('swc_q', 0.6, -1, -99.0,   0.038,   2.0, 0.3971,  1.88, 13, 100, 'mgamma', -99.0, -99.0,   1.0,    1.0,  5.13e-5, 2.294e-2, 'ss-rayleigh-gans', 'corPowerLaw_5.511054_0.25', -99.0),
       ('gwc_q', 1.0, -1, -99.0,  500.86,  3.18,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,  5.37,   1.06,  2.11e-4,   1.3e-2,       'mie-sphere',   'khvorostyanov01_spheres', -99.0), 
       ('hwc_q', 1.0, -1, -99.0,  392.33,   3.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   5.0,    1.0,  1.87e-4,   1.1e-2,       'mie-sphere',   'khvorostyanov01_spheres', -99.0)],
      dtype=[('hydro_name', 'S15'), ('as_ratio', '<f8'), ('liq_ice', '<i8'), ('rho_ms', '<f8'), ('a_ms', '<f8'), ('b_ms', '<f8'), ('alpha_as', '<f8'), ('beta_as', '<f8'), ('moment_in', '<i8'), ('nbin', '<i8'), ('dist_name', 'S15'), ('p_1', '<f8'), ('p_2', '<f8'), ('p_3', '<f8'), ('p_4', '<f8'), ('d_1', '<f8'), ('d_2', '<f8'), ('scat_name', 'S20'), ('vel_size_mod', 'S30'), ('canting', '<f8')]
      )

#def runHatpro(datestr):
# Meteogram
ICON_filename = ICON_folder + 'METEOGRAM_patch001_' + datestr + '_joyce.nc'
pam = pyPamtra.importer.readIcon2momMeteogram(ICON_filename,
											  descriptorFile,
											  timeidx=np.arange(990,991,1),
											  verbosity=1)


QC = 1.0*pam.p['hydro_q'][...,0]
QI = 1.0*pam.p['hydro_q'][...,1]
QR = 1.0*pam.p['hydro_q'][...,2]
QS = 1.0*pam.p['hydro_q'][...,3]
QG = 1.0*pam.p['hydro_q'][...,4]
QH = 1.0*pam.p['hydro_q'][...,5]

QNC = 1.0*pam.p['hydro_n'][...,0]
QNI = 1.0*pam.p['hydro_n'][...,1]
QNR = 1.0*pam.p['hydro_n'][...,2]
QNS = 1.0*pam.p['hydro_n'][...,3]
QNG = 1.0*pam.p['hydro_n'][...,4]
QNH = 1.0*pam.p['hydro_n'][...,5]


# SETTINGS
pam.nmlSet['active'] = True
pam.nmlSet["radar_mode"] = "moments"
pam.nmlSet['passive'] = False # Passive is time consuming
pam.set['verbose'] = 0 # set verbosity levels
pam.set['pyVerbose'] = 1 # change to 0 if you do not want to see job progress number
pam.p['turb_edr'][:] = 1.0e-4
pam.nmlSet['radar_airmotion'] = True
pam.nmlSet['radar_airmotion_vmax'] = 10.0
pam.nmlSet['radar_airmotion_vmin'] = -10.0

# X-band specific parameters Joyrad10
frequency = 9.6
pam.nmlSet['radar_fwhr_beamwidth_deg'] = 1.0
pam.nmlSet['radar_integration_time'] = 2.0
pam.nmlSet['radar_k2'] = 0.93
pam.nmlSet['radar_max_v'] = 78.07291
pam.nmlSet['radar_min_v'] = -78.07291
pam.nmlSet['radar_nfft'] = 4096
pam.nmlSet['radar_no_ave'] = 10
pam.nmlSet['radar_pnoise0'] = -48.0

pam.runParallelPamtra(np.array([frequency]), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF('/data/optimice/pamtra_runs/tripex-pol/data/'+datestr+'all_hydro_mom_X.nc')

pam.p['hydro_n'][...,3] = 0.0
pam.p['hydro_q'][...,3] = 0.0
pam.runParallelPamtra(np.array([frequency]), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF('/data/optimice/pamtra_runs/tripex-pol/data/'+datestr+'no_snow_mom_X.nc')

pam.p['hydro_n'][:] = 0.0
pam.p['hydro_q'][:] = 0.0
pam.p['hydro_n'][...,3] = 1.0*QNS
pam.p['hydro_q'][...,3] = 1.0*QS
pam.runParallelPamtra(np.array([frequency]), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF('/data/optimice/pamtra_runs/tripex-pol/data/'+datestr+'only_snow_mom_X.nc')

pam.p['hydro_n'][:] = 0.0
pam.p['hydro_q'][:] = 0.0
pam.p['hydro_n'][...,1] = 1.0*QNI
pam.p['hydro_q'][...,1] = 1.0*QI
pam.runParallelPamtra(np.array([frequency]), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF('/data/optimice/pamtra_runs/tripex-pol/data/'+datestr+'only_ice_mom_X.nc')

pam.p['hydro_n'][:] = 0.0
pam.p['hydro_q'][:] = 0.0
pam.p['hydro_n'][...,0] = 1.0*QNC
pam.p['hydro_q'][...,0] = 1.0*QC
pam.p['hydro_n'][...,2] = 1.0*QNR
pam.p['hydro_q'][...,2] = 1.0*QR
pam.runParallelPamtra(np.array([frequency]), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF('/data/optimice/pamtra_runs/tripex-pol/data/'+datestr+'only_liquid_mom_X.nc')

pam.p['hydro_n'][:] = 0.0
pam.p['hydro_q'][:] = 0.0
pam.p['hydro_n'][...,4] = 1.0*QNG
pam.p['hydro_q'][...,4] = 1.0*QG
pam.p['hydro_n'][...,5] = 1.0*QNH
pam.p['hydro_q'][...,5] = 1.0*QH
pam.runParallelPamtra(np.array([frequency]), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF('/data/optimice/pamtra_runs/tripex-pol/data/'+datestr+'only_graupel_hail_mom_X.nc')

# Ka-band specific parameters Joyrad35
frequency = 35.6
pam.nmlSet['radar_fwhr_beamwidth_deg'] = 0.6
pam.nmlSet['radar_integration_time'] = 2.0
pam.nmlSet['radar_k2'] = 0.93
pam.nmlSet['radar_max_v'] = 10.56824
pam.nmlSet['radar_min_v'] = -10.56824
pam.nmlSet['radar_nfft'] = 512
pam.nmlSet['radar_no_ave'] = 38
pam.nmlSet['radar_pnoise0'] = -64.0

pam.runParallelPamtra(np.array([frequency]), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
	#pam.writeResultsToNetCDF('/data/optimice/pamtra_runs/tripex-pol/data/'+datestr+'hatpro.nc') # SAVE OUTPUT

# W-band specific parameters Joyrad94
frequency = 94.0
pam.nmlSet['radar_fwhr_beamwidth_deg'] = 0.5
pam.nmlSet['radar_integration_time'] = 1.0
pam.nmlSet['radar_k2'] = 0.74
pam.nmlSet['radar_max_v'] = 6.8
pam.nmlSet['radar_min_v'] = -6.8
pam.nmlSet['radar_nfft'] = 512
pam.nmlSet['radar_no_ave'] = 17
pam.nmlSet['radar_pnoise0'] = -54.0

pam.runParallelPamtra(np.array([frequency]), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
	#pam.writeResultsToNetCDF('/data/optimice/pamtra_runs/tripex-pol/data/'+datestr+'hatpro.nc') # SAVE OUTPUT

# if (args.date[0] == 'all'):
#     alldates = [os.path.basename(x)[:8] for x in glob('/data/optimice/pamtra_runs/tripex-pol/data/*all_hydro_mom.nc')]
#     print('Running over all dates for which we have radar simulations',alldates)
#     for dd in alldates:
#         runHatpro(dd)
# else:
#     runHatpro(datestr = args.date[0])
