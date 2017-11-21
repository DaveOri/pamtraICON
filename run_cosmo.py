from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pn
import pyPamtra



descriptorFile = np.array([
      #['hydro_name' 'as_ratio' 'liq_ice' 'rho_ms' 'a_ms' 'b_ms' 'alpha_as' 'beta_as' 'moment_in' 'nbin' 'dist_name' 'p_1' 'p_2' 'p_3' 'p_4' 'd_1' 'd_2' 'scat_name' 'vel_size_mod' 'canting']
       ('cwc_q', -99.0, 1, -99.0, -99.0, -99.0, -99.0, -99.0, 3, 1, 'mono', -99.0, -99.0, -99.0, -99.0, 2e-05, -99.0, 'mie-sphere', 'khvorostyanov01_drops', -99.0),
       ('iwc_q', -99.0, -1, -99.0, 130.0, 3.0, 0.684, 2.0, 3, 1, 'mono_cosmo_ice', -99.0, -99.0, -99.0, -99.0, -99.0, -99.0, 'mie-sphere', 'heymsfield10_particles', -99.0),
       ('rwc_q', -99.0, 1, -99.0, -99.0, -99.0, -99.0, -99.0, 3, 50, 'exp', -99.0, -99.0, 8000000.0, -99.0, 0.00012, 0.006, 'mie-sphere', 'khvorostyanov01_drops', -99.0),
       ('swc_q', -99.0, -1, -99.0, 0.038, 2.0, 0.3971, 1.88, 3, 50, 'exp_cosmo_snow', -99.0, -99.0, -99.0, -99.0, 5.1e-11, 0.02, 'mie-sphere', 'heymsfield10_particles', -99.0),
       ('gwc_q', -99.0, -1, -99.0, 169.6, 3.1, -99.0, -99.0, 3, 50, 'exp', -99.0, -99.0, 4000000.0, -99.0, 1e-10, 0.01, 'mie-sphere', 'khvorostyanov01_spheres', -99.0)], 
      dtype=[('hydro_name', 'S15'), ('as_ratio', '<f8'), ('liq_ice', '<i8'), ('rho_ms', '<f8'), ('a_ms', '<f8'), ('b_ms', '<f8'), ('alpha_as', '<f8'), ('beta_as', '<f8'), ('moment_in', '<i8'), ('nbin', '<i8'), ('dist_name', 'S15'), ('p_1', '<f8'), ('p_2', '<f8'), ('p_3', '<f8'), ('p_4', '<f8'), ('d_1', '<f8'), ('d_2', '<f8'), ('scat_name', 'S15'), ('vel_size_mod', 'S30'), ('canting', '<f8')]
      )

constantFields = "/home/mech/workspace/pamtra/doc/tutorials/data/cosmo_constant_fields.nc"
fname='/home/mech/workspace/pamtra/doc/tutorials/data/LMK_gop9_test_fields_SynSatMic_201307241200-0900.nc'
pam = pyPamtra.importer.readCosmoDe1MomDataset(fname,"gop_fields_SynSatMic",descriptorFile,colIndex=0,verbosity=1,constantFields=constantFields)

#turn off passive calculations
pam.nmlSet["passive"] = False
#used frequencies
frequencies = [94]

#show some messages
pam.set["verbose"] = 0
pam.set["pyVerbose"] = 0

filterPam = np.zeros(pam._shape2D,dtype=bool)
filterPam[3:400,240] = True
print pam._shape3D
pam.filterProfiles(filterPam)
print pam._shape3D

pam.nmlSet['radar_attenuation'] = 'top-down'

pam.runParallelPamtra(frequencies, pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=1)
#pam.runPamtra([94])

outName = "cosmo_field_20130724_94.nc"
pam.writeResultsToNetCDF(outName, ncForm='NETCDF3_CLASSIC')

frequency = 0
radar_npol = 0
gridy = 0
radar_npeak = 0
Ze = pam.r["Ze"][:,gridy,:,frequency,radar_npol,radar_npeak]
H = pam.p["hgt"][:,gridy]
lon = pam.p["lon"][:,gridy]

matplotlib.rcParams['figure.figsize'] = (10.0, 8.0)
plt.pcolormesh(lon,H.T,Ze.T,vmin=-50,vmax=30)
plt.ylim(0,np.max(H))
plt.xlim(np.min(lon),np.max(lon))
plt.xlabel("longitude [deg]")
plt.ylabel("height [m]")
plt.colorbar(label="Ze [dBz]")
plt.show()
