import numpy as np
import pyPamtra

# Initialize PyPAMTRA instance
pam = pyPamtra.pyPamtra()

print "Settings"
descriptorFile = np.array([
      # ['hydro_name' 'as_ratio' 'liq_ice' 'rho_ms' 'a_ms' 'b_ms' 'alpha_as' 'beta_as' 'moment_in' 'nbin' 'dist_name' 'p_1' 'p_2' 'p_3' 'p_4' 'd_1' 'd_2' 'scat_name' 'vel_size_mod' 'canting']
       ('cwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,   2.0e-6,   8.0e-5, 'mie-sphere', 'corAtlas_9.292_9.623_622.2', -99.0),
       ('iwc_q', 1.0, -1, -99.0, 1.58783,  2.56,  0.684,   2.0, 13, 100, 'mgamma', -99.0, -99.0, 1.564, 0.8547, 1.744e-5, 9.369e-3, 'ssrg-rt3_0.18_0.89_2.06_0.08',   'corPowerLaw_30.606_0.5533',  -99.0),
       ('rwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,  0.00012,   8.2e-3, 'mie-sphere', 'corAtlas_9.292_9.623_622.2',  -99.0),
       ('swc_q', 0.6, -1, -99.0,   0.038,   2.0, 0.3971,  1.88, 13, 100, 'mgamma', -99.0, -99.0,   1.0,    1.0,  5.13e-5, 2.294e-2, 'ssrg-rt3_0.25_1.00_1.66_0.04',   'corPowerLaw_5.511054_0.25',  -99.0),
       ('gwc_q', 1.0, -1, -99.0,  500.86,  3.18,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,  5.37,   1.06,  2.11e-4,   1.3e-2, 'mie-sphere', 'corPowerLaw_406.67_0.85',    -99.0), 
       ('hwc_q', 1.0, -1, -99.0,  392.33,   3.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   5.0,    1.0,  1.87e-4,   1.1e-2, 'mie-sphere', 'corPowerLaw_106.33_0.5',     -99.0)],
      dtype=[('hydro_name', 'S15'), ('as_ratio', '<f8'), ('liq_ice', '<i8'), ('rho_ms', '<f8'), ('a_ms', '<f8'), ('b_ms', '<f8'), ('alpha_as', '<f8'), ('beta_as', '<f8'), ('moment_in', '<i8'), ('nbin', '<i8'), ('dist_name', 'S15'), ('p_1', '<f8'), ('p_2', '<f8'), ('p_3', '<f8'), ('p_4', '<f8'), ('d_1', '<f8'), ('d_2', '<f8'), ('scat_name', 'S30'), ('vel_size_mod', 'S30'), ('canting', '<f8')]
      )
#descriptorFile = np.array([
#      # ['hydro_name' 'as_ratio' 'liq_ice' 'rho_ms' 'a_ms' 'b_ms' 'alpha_as' 'beta_as' 'moment_in' 'nbin' 'dist_name' 'p_1' 'p_2' 'p_3' 'p_4' 'd_1' 'd_2' 'scat_name' 'vel_size_mod' 'canting']
#       ('cwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,   2.0e-6,   8.0e-5, 'mie-sphere', 'khvorostyanov01_drops', -99.0),
#       ('iwc_q', 1.0, -1, -99.0, 1.58783,  2.56,  0.684,   2.0, 13, 100, 'mgamma', -99.0, -99.0, 1.564, 0.8547, 1.744e-5, 9.369e-3, 'ssrg-rt3',   'powerLaw_30.606_0.5533',  -99.0),
#       ('rwc_q', 1.0, -1, 916.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,  0.00012,   8.2e-3, 'mie-sphere', 'powerLaw_494.74_0.7031',  -99.0),
#       ('swc_q', 0.6, -1, -99.0,   0.038,   2.0, 0.3971,  1.88, 13, 100, 'mgamma', -99.0, -99.0,   1.0,    1.0,  5.13e-5, 2.294e-2, 'ssrg-rt3',   'powerLaw_5.511054_0.25',  -99.0),
#       ('gwc_q', 1.0, -1, -99.0,  500.86,  3.18,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,  5.37,   1.06,  2.11e-4,   1.3e-2, 'mie-sphere', 'powerLaw_406.67_0.85',    -99.0), 
#       ('hwc_q', 1.0, -1, -99.0,  392.33,   3.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   5.0,    1.0,  1.87e-4,   1.1e-2, 'mie-sphere', 'powerLaw_106.33_0.5',     -99.0)],
#      dtype=[('hydro_name', 'S15'), ('as_ratio', '<f8'), ('liq_ice', '<i8'), ('rho_ms', '<f8'), ('a_ms', '<f8'), ('b_ms', '<f8'), ('alpha_as', '<f8'), ('beta_as', '<f8'), ('moment_in', '<i8'), ('nbin', '<i8'), ('dist_name', 'S15'), ('p_1', '<f8'), ('p_2', '<f8'), ('p_3', '<f8'), ('p_4', '<f8'), ('d_1', '<f8'), ('d_2', '<f8'), ('scat_name', 'S20'), ('vel_size_mod', 'S30'), ('canting', '<f8')]
#      )

#descriptorFile = np.array([
#      #['hydro_name' 'as_ratio' 'liq_ice' 'rho_ms' 'a_ms' 'b_ms' 'alpha_as' 'beta_as' 'moment_in' 'nbin' 'dist_name' 'p_1' 'p_2' 'p_3' 'p_4' 'd_1' 'd_2' 'scat_name' 'vel_size_mod' 'canting']
#       ('cwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,   2.0e-6,   8.0e-5,       'mie-sphere', 'khvorostyanov01_drops', -99.0),
#       ('iwc_q', 1.0, -1, -99.0, 1.58783,  2.56,  0.684,   2.0, 13, 100, 'mgamma', -99.0, -99.0, 1.564, 0.8547, 1.744e-5, 9.369e-3, 'ss-rayleigh-gans', 'corPowerLaw_30.606_0.5533', -99.0),
#       ('rwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,  0.00012,   8.2e-3,       'mie-sphere', 'khvorostyanov01_drops', -99.0),
#       ('swc_q', 0.6, -1, -99.0,   0.038,   2.0, 0.3971,  1.88, 13, 100, 'mgamma', -99.0, -99.0,   1.0,    1.0,  5.13e-5, 2.294e-2, 'ss-rayleigh-gans', 'corPowerLaw_5.511054_0.25', -99.0),
#       ('gwc_q', 1.0, -1, -99.0,  500.86,  3.18,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,  5.37,   1.06,  2.11e-4,   1.3e-2,       'mie-sphere', 'khvorostyanov01_spheres', -99.0), 
#       ('hwc_q', 1.0, -1, -99.0,  392.33,   3.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   5.0,    1.0,  1.87e-4,   1.1e-2,       'mie-sphere', 'khvorostyanov01_spheres', -99.0)],
#      dtype=[('hydro_name', 'S15'), ('as_ratio', '<f8'), ('liq_ice', '<i8'), ('rho_ms', '<f8'), ('a_ms', '<f8'), ('b_ms', '<f8'), ('alpha_as', '<f8'), ('beta_as', '<f8'), ('moment_in', '<i8'), ('nbin', '<i8'), ('dist_name', 'S15'), ('p_1', '<f8'), ('p_2', '<f8'), ('p_3', '<f8'), ('p_4', '<f8'), ('d_1', '<f8'), ('d_2', '<f8'), ('scat_name', 'S20'), ('vel_size_mod', 'S30'), ('canting', '<f8')]
#      )
for hyd in descriptorFile:
	pam.df.addHydrometeor(hyd)

hydrometcount = len(descriptorFile)
# def mixing ratio and number concentration to simulate a wide part of q/N
N = 100
mix_rat = np.logspace(-8, -4, N)
num_conc = np.ones(mix_rat.shape)*1.0e2  # for the DWRs just q/N matters

pam.nmlSet["passive"] = False  # Activate this for Microwave radiometer
pam.nmlSet["active"] = True    # Activate this for Cloud radar
pam.nmlSet["radar_mode"] = 'spectrum'
pam.nmlSet["liq_mod"] = 'TKC'
# pam.nmlSet["save_psd"] = False    # save particle size distribution
# pam.nmlSet["radar_attenuation"] = "disabled" #"bottom-up"
# pam.nmlSet["hydro_adaptive_grid"] = False    # uncomment (or set to false) this line before changing max and min diameter of size distribution
# pam.nmlSet["conserve_mass_rescale_dsd"] = False    #necessary because of separating into m-D-relationship regions
# show some messages
# pam.nmlSet["hydro_fullspec"] = True #use full-spectra as input
# pam.nmlSet["radar_allow_negative_dD_dU"] = True #allow negative dU dD which can happen at the threshold between different particle species
pam.nmlSet["radar_nfft"] = 2048
pam.nmlSet["radar_no_ave"] = 60
pam.nmlSet["radar_pnoise0"] = -164.0
pam.nmlSet["radar_max_v"] = 10.0
pam.nmlSet["radar_min_v"] = -10.0
pam.nmlSet["radar_integration_time"] = 5.0
pam.set["verbose"] = 0
pam.set["pyVerbose"] = 0
# Generate PAMTRA data dictonary
pamData = dict()
vec_shape = [hydrometcount, N]
vec_shape_hyd = [hydrometcount, N, hydrometcount]
# define athmospheric variables
pamData["press"] = 73500*np.ones(vec_shape)  # Pressure
pamData["relhum"] = 90*np.ones(vec_shape)  # Relative Humidity
# pamData["timestamp"] = 1479945600*np.ones(vec_shape) #iconData["time"]
# Juelich coordinates
# pamData["lat"] = 50.901*np.ones(vec_shape[0])
# pamData["lon"] = 6.391*np.ones(vec_shape[0])
pamData["lfrac"] = 1.0*np.ones(vec_shape[0])  #icon_frland[timeframe]
# pamData["wind10u"] = 1*np.ones(vec_shape[0]) #iconData["u"][0,1]
# pamData["wind10v"] = 1*np.ones(vec_shape[0]) #iconData["v"][0,1]
pamData["hgt"] = np.ones(vec_shape)*np.arange(10.,1200.,(1200.-0.)/N)
pamData["temp"] = 268*np.ones(vec_shape)  # temperature
pamData["hydro_q"] = np.zeros(vec_shape_hyd)  # TODO: query number of categories (here 6) automatically
pamData["hydro_n"] = np.zeros(vec_shape_hyd)  # np.zeros([1,num_lev,4]) #TODO: query number of categories (here 4) automatically

print "Creating profile"
for hyd in range(hydrometcount):
    pamData["hydro_q"][hyd, :, hyd] = 1.0*mix_rat
    pamData["hydro_n"][hyd, :, hyd] = 1.0*num_conc
# Add them to pamtra object and create profile
pam.createProfile(**pamData)
# set hydrometeor properties

pam.runPamtra([9.4, 35.5, 94.])
pam.writeResultsToNetCDF('data/idealized_hydro.nc')

#pam.nmlSet["radar_mode"] = 'moments'
#pam.runPamtra([9.4, 35.5, 94.])
#pam.writeResultsToNetCDF('data/idealized_hydro_spectrum.nc')

pam.nmlSet["radar_mode"] = 'simple'
pam.set["verbose"] = 0
pam.runPamtra([9.4, 35.5, 94.])
pam.writeResultsToNetCDF('data/idealized_hydro_simple.nc')
