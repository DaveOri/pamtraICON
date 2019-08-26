from __future__ import division
import pyPamtra
import numpy as np
import argparse
from radar_settings import radarlib, hydrodict

parser =  argparse.ArgumentParser(description='do plots for QuickLookBrowser')
parser.add_argument('-d','--date', nargs=1,
	                help='gimme datestring in the format YYYYMMDD')
parser.add_argument('-r','--radarset', nargs=1,
	                help='gimme radar name',
	                choices=radarlib.keys())
parser.add_argument('-hy','--hydroset', nargs=1,
	                help='gimme hydrosettings',
	                choices=hydrodict.keys())
parser.add_argument('-dt','--datapath', nargs=1,
	                help='gimme fulldatapath for saving output')
parser.add_argument('-p', '--patch', nargs=1, help='tell me the 3 padded patch number like 001')
parser.add_argument('-np', '--numproc', nargs=1, help='gimme the number of processors I can use default = 8')
parser.print_help()
args = parser.parse_args()
datestr = args.date[0]
radarstr = args.radarset[0]
hydrostr = args.hydroset[0]
print datestr, radarstr, hydrostr

cores = 8 # number of parallel cores
if args.numproc is not None:
	cores = int(args.numproc[0])

# PATHS
datapath = '/data/optimice/pamtra_runs/tripex-pol/data/'
if args.datapath is not None:
	datapath = args.datapath[0]
# Directory where all the .nc meteograms are stored
ICON_folder = '/data/inscape/icon/experiments/juelich/testbed/testbed_' + datestr + '/'

mom1 = False
if '_1mom' in datestr:
	print datestr
	mom1 = True
datestr = datestr[:8]

patch = '001'
patchsave = False
if args.patch is not None:
	patch = args.patch[0]
	patchsave = True

# Folder where your descriptor files are stored (you can use mine, Mario's or the default pamtra)
# descriptor_folder = '/home/dori/descriptorfiles/'

# INIT
#descriptor_filename = 'descriptor_file_2m_ssrgNEWpowerLaw.txt'
#descriptorFile2mom = np.array([
#      #['hydro_name' 'as_ratio' 'liq_ice' 'rho_ms' 'a_ms' 'b_ms' 'alpha_as' 'beta_as' 'moment_in' 'nbin' 'dist_name' 'p_1' 'p_2' 'p_3' 'p_4' 'd_1' 'd_2' 'scat_name' 'vel_size_mod' 'canting']
#       ('cwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,   2.0e-6,   8.0e-5, 'mie-sphere', 'corPowerLaw_24388657.6_2.0', -99.0),
#       ('iwc_q', 1.0, -1, -99.0, 1.58783,  2.56,  0.684,   2.0, 13, 100, 'mgamma', -99.0, -99.0, 1.564, 0.8547, 1.744e-5, 9.369e-3, 'ssrg-rt3',   'corPowerLaw_30.606_0.5533',  -99.0),
#       ('rwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,  0.00012,   8.2e-3, 'mie-sphere', 'corPowerLaw_494.74_0.7031',  -99.0),
#       ('swc_q', 0.6, -1, -99.0,   0.038,   2.0, 0.3971,  1.88, 13, 100, 'mgamma', -99.0, -99.0,   1.0,    1.0,  5.13e-5, 2.294e-2, 'ssrg-rt3',   'corPowerLaw_5.511054_0.25',  -99.0),
#       ('gwc_q', 1.0, -1, -99.0,  500.86,  3.18,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,  5.37,   1.06,  2.11e-4,   1.3e-2, 'mie-sphere', 'corPowerLaw_406.67_0.85',    -99.0), 
#       ('hwc_q', 1.0, -1, -99.0,  392.33,   3.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   5.0,    1.0,  1.87e-4,   1.1e-2, 'mie-sphere', 'corPowerLaw_106.33_0.5',     -99.0)],
#      dtype=[('hydro_name', 'S15'), ('as_ratio', '<f8'), ('liq_ice', '<i8'), ('rho_ms', '<f8'), ('a_ms', '<f8'), ('b_ms', '<f8'), ('alpha_as', '<f8'), ('beta_as', '<f8'), ('moment_in', '<i8'), ('nbin', '<i8'), ('dist_name', 'S15'), ('p_1', '<f8'), ('p_2', '<f8'), ('p_3', '<f8'), ('p_4', '<f8'), ('d_1', '<f8'), ('d_2', '<f8'), ('scat_name', 'S20'), ('vel_size_mod', 'S30'), ('canting', '<f8')]
#      )
descriptorFile2mom = np.array([
      # ['hydro_name' 'as_ratio' 'liq_ice' 'rho_ms' 'a_ms' 'b_ms' 'alpha_as' 'beta_as' 'moment_in' 'nbin' 'dist_name' 'p_1' 'p_2' 'p_3' 'p_4' 'd_1' 'd_2' 'scat_name' 'vel_size_mod' 'canting']
       ('cwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,   2.0e-6,   8.0e-5, 'mie-sphere',                     'corAtlas_9.292_9.623_622.2', -99.0),
       ('iwc_q', 1.0, -1, -99.0, 1.58783,  2.56,  0.684,   2.0, 13, 100, 'mgamma', -99.0, -99.0, 1.564, 0.8547, 1.744e-5, 9.369e-3, 'ssrg-rt3_0.18_0.89_2.06_0.08',   'corPowerLaw_30.606_0.5533',  -99.0),
       ('rwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,  0.00012,   8.2e-3, 'mie-sphere',                     'corAtlas_9.292_9.623_622.2',  -99.0),
       ('swc_q', 0.6, -1, -99.0,   0.038,   2.0, 0.3971,  1.88, 13, 100, 'mgamma', -99.0, -99.0,   1.0,    1.0,  5.13e-5, 2.294e-2, 'ssrg-rt3_0.25_1.00_1.66_0.04',   'corPowerLaw_5.511054_0.25',  -99.0),
       ('gwc_q', 1.0, -1, -99.0,  500.86,  3.18,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,  5.37,   1.06,  2.11e-4,   1.3e-2, 'mie-sphere',                     'corPowerLaw_406.67_0.85',    -99.0), 
       ('hwc_q', 1.0, -1, -99.0,  392.33,   3.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   5.0,    1.0,  1.87e-4,   1.1e-2, 'mie-sphere',                     'corPowerLaw_106.33_0.5',     -99.0)],
      dtype=[('hydro_name', 'S15'), ('as_ratio', '<f8'), ('liq_ice', '<i8'), ('rho_ms', '<f8'), ('a_ms', '<f8'), ('b_ms', '<f8'), ('alpha_as', '<f8'), ('beta_as', '<f8'), ('moment_in', '<i8'), ('nbin', '<i8'), ('dist_name', 'S15'), ('p_1', '<f8'), ('p_2', '<f8'), ('p_3', '<f8'), ('p_4', '<f8'), ('d_1', '<f8'), ('d_2', '<f8'), ('scat_name', 'S30'), ('vel_size_mod', 'S30'), ('canting', '<f8')]
      )

descriptorFile1mom = np.array([ # TODO to be reviewed, coefficients for m-D and v-D are changing
  #['hydro_name' 'as_ratio' 'liq_ice' 'rho_ms' 'a_ms' 'b_ms' 'alpha_as' 'beta_as' 'moment_in' 'nbin' 'dist_name' 'p_1' 'p_2' 'p_3' 'p_4' 'd_1' 'd_2' 'scat_name' 'vel_size_mod' 'canting']
   ('cwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0,  3,  1,   'mono',           -99.0, -99.0, -99.0, -99.0,  2.0e-5,  -99.0, 'mie-sphere', 'corPowerLaw_24388657.6_2.0', -99.0),
   ('iwc_q', 0.2, -1, -99.0,   130.0,   3.0,  0.684,   2.0,  3,  1,   'mono_cosmo_ice', -99.0, -99.0, -99.0, -99.0,   -99.0,  -99.0, 'ssrg-rt3_0.18_0.89_2.06_0.08', 'corPowerLaw_30.606_0.5533', -99.0),
   ('rwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0,  3,  100, 'exp',            -99.0, -99.0, 8.0e6, -99.0,  1.2e-4, 6.0e-3, 'mie-sphere', 'corPowerLaw_130.0_0.5', -99.0),
   ('swc_q', 0.6, -1, -99.0,   0.038,   2.0, 0.3971,  1.88,  3,  100, 'exp_cosmo_snow', -99.0, -99.0, -99.0, -99.0, 5.1e-11, 1.0e-2, 'ssrg-rt3_0.25_1.00_1.66_0.04', 'corPowerLaw_4.9_0.25', -99.0),
   ('gwc_q', 1.0, -1, -99.0,   169.6,   3.1,  -99.0, -99.0,  3,  100, 'exp',            -99.0, -99.0, 4.0e6, -99.0, 1.0e-10, 1.0e-2, 'mie-sphere', 'corPowerLaw_406.67_0.85', -99.0)],
  dtype=[('hydro_name', 'S15'), ('as_ratio', '<f8'), ('liq_ice', '<i8'), ('rho_ms', '<f8'), ('a_ms', '<f8'), ('b_ms', '<f8'), ('alpha_as', '<f8'), ('beta_as', '<f8'), ('moment_in', '<i8'), ('nbin', '<i8'), ('dist_name', 'S15'), ('p_1', '<f8'), ('p_2', '<f8'), ('p_3', '<f8'), ('p_4', '<f8'), ('d_1', '<f8'), ('d_2', '<f8'), ('scat_name', 'S20'), ('vel_size_mod', 'S30'), ('canting', '<f8')]
  )

# Meteogram
ICON_filename = ICON_folder + 'METEOGRAM_patch' + patch +'_' + datestr + '_joyce.nc'
print ICON_filename
if mom1:
	pam = pyPamtra.importer.readIcon1momMeteogram(ICON_filename,
	    										  descriptorFile1mom,
		    									  timeidx=None,#np.arange(1200,2400),
			    								  verbosity=1,
				    							  hydro_content=hydrodict[hydrostr][:5])
else:
	pam = pyPamtra.importer.readIcon2momMeteogram(ICON_filename,
					    						  descriptorFile2mom,
						    					  timeidx=None,#np.arange(1200,2400),
							    				  verbosity=1,
								    			  hydro_content=hydrodict[hydrostr])

# SETTINGS
pam.nmlSet['active'] = True
pam.nmlSet["radar_mode"] = 'moments'#"spectrum"
pam.nmlSet['passive'] = False # Passive is time consuming
pam.set['verbose'] = 0 # set verbosity levels
pam.set['pyVerbose'] = 1 # change to 0 if you do not want to see job progress number
pam.p['turb_edr'][:] = 1.0e-4
pam.nmlSet['radar_airmotion'] = True
pam.nmlSet['radar_airmotion_vmin'] = 0.0 # workaround to potential bug in radar_spectrum
pam.nmlSet['radar_airmotion_model'] = 'constant'

def set_radar_properties(pam,radarlib,radar):
	radarbook = radarlib[radar]
	for k in radarbook.keys():
		print k
		if 'radar' in k: # avoid to set frequency in the nmlSet
			print k, radarbook[k]
			pam.nmlSet[k] = radarbook[k]
	return pam, radarbook['frequency']

def run_radar_simulation(pam, radarname, hydroconf):
	pam, frequency = set_radar_properties(pam, radarlib, radarname)
	pam.runParallelPamtra(np.array([frequency]), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
	if patchsave:
		pam.writeResultsToNetCDF(datapath+hydroconf+'/'+datestr+hydroconf+patch+'_'+pam.nmlSet["radar_mode"][:3]+'_'+radarname+'.nc')
	else:
		pam.writeResultsToNetCDF(datapath+hydroconf+'/'+datestr+hydroconf+'_'+pam.nmlSet["radar_mode"][:3]+'_'+radarname+'.nc')
	return pam

def runHatpro(datestr):
	# Meteogram
	ICON_filename = ICON_folder + 'METEOGRAM_patch001_' + datestr + '_joyce.nc'
	pam = pyPamtra.importer.readIcon2momMeteogram(ICON_filename,
		                                          descriptorFile2mom,
		                                          timeidx=None,
		                                          verbosity=1,
		                                          hydro_content=hydrodict[hydrostr])
	# SETTINGS
	pam.nmlSet['active'] = False
	pam.nmlSet['passive'] = True # Passive is time consuming
	pam.set['verbose'] = 0 # set verbosity levels
	pam.set['pyVerbose'] = 1 # change to 0 if you do not want to see job progress number

	f_hatpro_Kband = [22.24, 23.04, 23.84, 25.44, 26.24, 27.84, 31.40]
	f_hatpro_Vband = [51.26, 52.28, 53.86, 54.94, 56.66, 57.30, 58.00]
	frequencies = f_hatpro_Kband + f_hatpro_Vband
	pam.runParallelPamtra(np.array(frequencies), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
	pam.writeResultsToNetCDF(datapath+datestr+'hatpro.nc') # SAVE OUTPUT

# RUN PAMTRA
if (radarstr == 'hatpro'):
	runHatpro(datestr=datestr)
else:
	pam = run_radar_simulation(pam, radarstr, hydrostr)
