from __future__ import division
import pyPamtra
import numpy as np
import argparse
from radar_settings import radarlib, hydrodict

cores = 8 # number of parallel cores
patches = ['001','002','003','004']

# INIT
descriptorFile = np.array([
      #['hydro_name' 'as_ratio' 'liq_ice' 'rho_ms' 'a_ms' 'b_ms' 'alpha_as' 'beta_as' 'moment_in' 'nbin' 'dist_name' 'p_1' 'p_2' 'p_3' 'p_4' 'd_1' 'd_2' 'scat_name' 'vel_size_mod' 'canting']
       ('cwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,   2.0e-6,   8.0e-5,       'mie-sphere', 'corPowerLaw_24388657.6_2.0', -99.0),
       ('iwc_q', 1.0, -1, -99.0, 1.58783,  2.56,  0.684,   2.0, 13, 100, 'mgamma', -99.0, -99.0, 1.564, 0.8547, 1.744e-5, 9.369e-3, 'ss-rayleigh-gans', 'corPowerLaw_30.606_0.5533',  -99.0),
       ('rwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,  0.00012,   8.2e-3,       'mie-sphere', 'corPowerLaw_494.74_0.7031',  -99.0),
       ('swc_q', 0.6, -1, -99.0,   0.038,   2.0, 0.3971,  1.88, 13, 100, 'mgamma', -99.0, -99.0,   1.0,    1.0,  5.13e-5, 2.294e-2, 'ss-rayleigh-gans', 'corPowerLaw_5.511054_0.25',  -99.0),
       ('gwc_q', 1.0, -1, -99.0,  500.86,  3.18,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,  5.37,   1.06,  2.11e-4,   1.3e-2,       'mie-sphere', 'corPowerLaw_406.67_0.85',    -99.0), 
       ('hwc_q', 1.0, -1, -99.0,  392.33,   3.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   5.0,    1.0,  1.87e-4,   1.1e-2,       'mie-sphere', 'corPowerLaw_106.33_0.5',     -99.0)],
      dtype=[('hydro_name', 'S15'), ('as_ratio', '<f8'), ('liq_ice', '<i8'), ('rho_ms', '<f8'), ('a_ms', '<f8'), ('b_ms', '<f8'), ('alpha_as', '<f8'), ('beta_as', '<f8'), ('moment_in', '<i8'), ('nbin', '<i8'), ('dist_name', 'S15'), ('p_1', '<f8'), ('p_2', '<f8'), ('p_3', '<f8'), ('p_4', '<f8'), ('d_1', '<f8'), ('d_2', '<f8'), ('scat_name', 'S20'), ('vel_size_mod', 'S40'), ('canting', '<f8')]
      )
# Meteogram
hydrostr = 'all_hydro'
for patch in patches:
  ICON_filename = '/data/inscape/icon/experiments/nyalesund/2016-11-05/METEOGRAM_patch'+patch+'_awipev.nc'
  pam = pyPamtra.importer.readIcon2momMeteogram(ICON_filename,
                                                descriptorFile,
                                                timeidx=None,#np.arange(1200,2400),
                                                verbosity=1,
                                                hydro_content=hydrodict[hydrostr])

  # SETTINGS
  pam.nmlSet['active'] = True
  pam.nmlSet["radar_mode"] = "spectrum"
  pam.nmlSet['passive'] = False # Passive is time consuming
  pam.set['verbose'] = 0 # set verbosity levels
  pam.set['pyVerbose'] = 1 # change to 0 if you do not want to see job progress number
  pam.p['turb_edr'][:] = 1.0e-4
  pam.nmlSet['radar_airmotion'] = True
  pam.nmlSet['radar_airmotion_vmin'] = 0.0 # workaround to potential bug in radar_spectrum
  pam.nmlSet['radar_airmotion_model'] = 'constant'

  frequency = [35.6,94]
  pam.runParallelPamtra(np.array(frequency), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
  pam.writeResultsToNetCDF('/data/optimice/pamtra_runs/nyalesund/20171105/METEOGRAM_awipev_patch'+patch+'.nc')