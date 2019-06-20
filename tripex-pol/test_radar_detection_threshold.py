# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 16:23:02 2019

@author: dori
"""

import pyPamtra
import matplotlib.pyplot as plt
import numpy as np
import copy

from radar_settings import radarlib
plt.close('all')

def hildebrand_sekhon(spectrum, nAvg):
  
  N = len(spectrum) - 1 # index of the last element
  sortSpec = np.sort(spectrum)
  sortSpec[~np.isfinite(sortSpec)] = 0.0
  sum2 = sortSpec[0:-1].cumsum()**2
  sumSquares= (sortSpec[0:-1]**2).cumsum()
  nPts=np.arange(N)+1.0
  sigma = nAvg*(nPts*sumSquares-sum2)
  a=int(max(nPts[sigma < sum2]))

  if(a < 10):
    peak_noise=max(spectrum)
  else:
    peak_noise=sortSpec[a-1]
  signal_detected=1.0*spectrum
  signal_detected[spectrum < peak_noise] = np.nan
  return peak_noise, signal_detected

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
filepath = '/data/inscape/icon/experiments/juelich/testbed/testbed_'
datestr = '20181124'
datestr = '20181219'
filename = filepath + datestr +'/METEOGRAM_patch001_'+datestr+'_joyce.nc'

pam = pyPamtra.importer.readIcon2momMeteogram(filename,
                                              descriptorFile=descriptorFile,
                                              timeidx=[4800], 
                                              hydro_content=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

def joseFit(x,a,b, c0=0.0):
    return c0 + a*np.log10(b*x)

coeff = {'Joyrad35':[19.785, 3.44e-7],
         'Joyrad10':[19.678, 2.55e-6],
         'KiXPol':[10.0, 2.19, -92.04],
         'Grarad94':[20.394, 1.294e-6] }


#pam.p['hydro_q'][:] = 0.0
#pam.p['hydro_n'][:] = 0.0
#pam.p['hydro_q'][0,0,:,2] = 1.0e-8
#pam.p['hydro_n'][0,0,:,2] = 1.0e2




radarname = 'Joyrad35'
#radarname = 'Joyrad10'
#radarname = 'Grarad94'
#radarname = 'KiXPol'

radar = radarlib[radarname]

for k in radar.keys():
    if 'radar' in k:
        pam.nmlSet[k] = radar[k]

#pam.nmlSet['radar_pnoise0'] = -48
#pam.nmlSet['radar_no_ave'] = 1000
#pam.nmlSet['radar_nfft'] = 512

# SETTINGS
pam.nmlSet['active'] = True
pam.nmlSet["radar_mode"] = "spectrum"
pam.nmlSet['passive'] = False # Passive is time consuming
pam.set['verbose'] = 0 # set verbosity levels
pam.set['pyVerbose'] = 0 # change to 0 if you do not want to see job progress number
pam.p['turb_edr'][:] = 1.0e-4
pam.nmlSet['radar_airmotion'] = True
#pam.nmlSet['radar_airmotion_vmin'] = 0.0 # workaround to potential bug in radar_spectrum
pam.nmlSet['radar_airmotion_model'] = 'constant'
#pam.nmlSet['radar_save_noise_corrected_spectra'] = True
#pam.nmlSet['radar_noise_distance_factor'] = 10.0
#pam.nmlSet['radar_peak_min_snr'] = 0.0
pam.nmlSet['radar_no_ave'] = 20

#pam.runPamtra(radar['frequency'])
pam2 = copy.deepcopy(pam)
pam.runParallelPamtra(radar['frequency'])

H = pam.p['hgt'][0,0,:]
Hr = pam.r['radar_hgt'][0,0,:]
Z = pam.r['Ze'][0,0,:,0,0,0]
Z[Z == -9999.] = np.nan
SNR = pam.r['radar_snr'][0,0,:,0,0,0]
plt.figure()
plt.plot(Z,H,lw=4,label='Ze moments')
pam.writeResultsToNetCDF('test_moments.nc')


pam2.nmlSet["radar_mode"] = "simple"
#pam.runPamtra(radar['frequency'])
pam2.runParallelPamtra(radar['frequency'])

Ze = pam2.r['Ze'][0,0,:,0,0,0]
Ze[Ze == -9999.] = np.nan
dopplervel = pam.r['radar_vel'][0,:]
##dv = 1.0#(dopplervel[1:]-dopplervel[0:-1])[0]
spectrogram = pam.r['radar_spectra'][0,0,:,0,0,:]
spectrogram[spectrogram == -9999.0] = np.nan
linspectrogram = 10.0**(0.1*spectrogram)
linspectrogram[~np.isfinite(linspectrogram)] = 0.0
Ze_spec = 10.0*np.log10(linspectrogram.sum(axis=1))

plt.plot(Ze,H, marker='x', label='Ze simple')
pam2.writeResultsToNetCDF('test_simple.nc')

Noise = Ze - SNR
plt.plot(Noise,H,lw=4,label='Ze simple - SNR moments')

Ncomp = pam.nmlSet['radar_pnoise0'] + 20.0*np.log10(H*0.001)
Njose = joseFit(H,*coeff[radarname])
plt.plot(Ncomp,H,'--',label='p0 + dB(Hgt)')
plt.plot(Njose,H,'.', label="jose' fit ")
plt.plot(Ze_spec,H,label="Spectrum integral")
plt.ylim([0,12000])
plt.grid()
plt.xlabel('dB')
plt.ylabel('height [m]')
plt.legend()
plt.savefig(radarname + 'test_reflectivity_threshold.png')

print np.nanmax(Noise-Z)

plt.figure()
plt.pcolormesh(dopplervel,H,spectrogram, vmin=-70, vmax=-30, cmap='jet')
plt.ylim([0,12000])
plt.colorbar()
plt.savefig('spectrogram.png')

hgtidx = 64
plt.figure()
noiselvl = dopplervel*0.0
noiselvl[:] = pam.nmlSet['radar_pnoise0'] + 20.0*np.log10(H[hgtidx]*0.001) -10.0*np.log10(pam.nmlSet['radar_nfft'])
spectrum = spectrogram[hgtidx,:]
linspectrum = 10.0**(0.1*spectrum)
peaknoise, detectedSpectrum = hildebrand_sekhon(linspectrum, pam.nmlSet['radar_no_ave'])
plt.plot(dopplervel,spectrogram[hgtidx,:], label='pamtra spectrum')
plt.plot(dopplervel, noiselvl, label = 'noise level pamtra')
plt.plot(dopplervel, noiselvl*0.0+10.0*np.log10(peaknoise), label='peak Noise HS')
plt.plot(dopplervel, 10.0*np.log10(detectedSpectrum),marker='x', label='detected Spectrum HS')
plt.xlim([-10,10])
plt.title('Hgt= ' + str(H[hgtidx])[:8])
plt.xlabel('Doppler velocity [m/s]')
plt.ylabel('Spectral power [dB]')
ZeHS = 10.0*np.log10(detectedSpectrum[np.isfinite(detectedSpectrum)].sum())

plt.legend(title='ZeHS= '+str(ZeHS)[:6]+' ZeP= '+str(Ze[hgtidx])[:6])
plt.savefig(radarname + 'spectrum.png')

