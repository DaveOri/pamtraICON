#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 16:30:09 2019

@author: dori
"""

import matplotlib.pyplot as plt
import netCDF4
import pandas as pd
import numpy as np

plt.rcParams.update({'font.size':14})

data = netCDF4.Dataset('../data/idealized_hydro.nc')
data_simple = netCDF4.Dataset('../data/idealized_hydro_simple.nc')
datavar = data.variables
Zvar = data_simple.variables['Ze']
Ze = Zvar[:]

Zc = Ze[0,0,:,:,0,0]
Zi = Ze[1,0,:,:,0,0]
Zr = Ze[2,0,:,:,0,0]
Zs = Ze[3,0,:,:,0,0]
Zg = Ze[4,0,:,:,0,0]
Zh = Ze[5,0,:,:,0,0]

Zc[:,0] -= 0.5
Zr[:,0] -= 0.5

descriptorFile = np.array([
       ('cwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,   2.0e-6,   8.0e-5, 'mie-sphere', 'corAtlas_9.292_9.623_622.2', -99.0),
       ('iwc_q', 1.0, -1, -99.0, 1.58783,  2.56,  0.684,   2.0, 13, 100, 'mgamma', -99.0, -99.0, 1.564, 0.8547, 1.744e-5, 9.369e-3, 'ssrg-rt3_0.18_0.89_2.06_0.08',   'corPowerLaw_30.606_0.5533',  -99.0),
       ('rwc_q', 1.0,  1, -99.0,   -99.0, -99.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   2.0,    1.0,  0.00012,   8.2e-3, 'mie-sphere', 'corAtlas_9.292_9.623_622.2',  -99.0),
       ('swc_q', 0.6, -1, -99.0,   0.038,   2.0, 0.3971,  1.88, 13, 100, 'mgamma', -99.0, -99.0,   1.0,    1.0,  5.13e-5, 2.294e-2, 'ssrg-rt3_0.25_1.00_1.66_0.04',   'corPowerLaw_5.511054_0.25',  -99.0),
       ('gwc_q', 1.0, -1, -99.0,  500.86,  3.18,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,  5.37,   1.06,  2.11e-4,   1.3e-2, 'mie-sphere', 'corPowerLaw_406.67_0.85',    -99.0), 
       ('hwc_q', 1.0, -1, -99.0,  392.33,   3.0,  -99.0, -99.0, 13, 100, 'mgamma', -99.0, -99.0,   5.0,    1.0,  1.87e-4,   1.1e-2, 'mie-sphere', 'corPowerLaw_106.33_0.5',     -99.0)],
      dtype=[('hydro_name', 'S15'), ('as_ratio', '<f8'), ('liq_ice', '<i8'), ('rho_ms', '<f8'), ('a_ms', '<f8'), ('b_ms', '<f8'), ('alpha_as', '<f8'), ('beta_as', '<f8'), ('moment_in', '<i8'), ('nbin', '<i8'), ('dist_name', 'S15'), ('p_1', '<f8'), ('p_2', '<f8'), ('p_3', '<f8'), ('p_4', '<f8'), ('d_1', '<f8'), ('d_2', '<f8'), ('scat_name', 'S30'), ('vel_size_mod', 'S30'), ('canting', '<f8')]
      )
cols = ['hydro_name', 'as_ratio', 'liq_ice', 'rho_ms', 'a_ms', 'b_ms',
        'alpha_as', 'beta_as', 'moment_in', 'nbin', 'dist_name',
        'p_1', 'p_2', 'p_3', 'p_4', 'd_1', 'd_2',
        'scat_name', 'vel_size_mod', 'canting']
df = pd.DataFrame(descriptorFile, columns=cols)
hydrometcount = len(descriptorFile)

# def mixing ratio and number concentration to simulate a wide part of q/N
N = 100
mix_rat = np.logspace(-8, -4, N)
num_conc = np.ones(mix_rat.shape)*1.0e2  # for the DWRs just q/N matters

avgMass = mix_rat/num_conc

avgMassc = 1.0e3*np.cbrt(1.0e-3*6*avgMass/np.pi)
avgMassi = 1.0e3*(avgMass/1.58783)**(1.0/2.56)
avgMassr = 1.0e3*np.cbrt(1.0e-3*6*avgMass/np.pi)
avgMasss = 1.0e3*(avgMass/0.038)**(1.0/2.0)
avgMassg = 1.0e3*(avgMass/500.86)**(1.0/3.18)
avgMassh = 1.0e3*(avgMass/392.33)**(1.0/3.0)


Zc = Ze[0,0,:,:,0,0]
Zi = Ze[1,0,:,:,0,0]
Zr = Ze[2,0,:,:,0,0]
Zs = Ze[3,0,:,:,0,0]
Zg = Ze[4,0,:,:,0,0]
Zh = Ze[5,0,:,:,0,0]

DWRc = np.array([Zc[:,0]-Zc[:,1],Zc[:,1]-Zc[:,2]])
DWRi = np.array([Zi[:,0]-Zi[:,1],Zi[:,1]-Zi[:,2]])
DWRr = np.array([Zr[:,0]-Zr[:,1],Zr[:,1]-Zr[:,2]])
DWRs = np.array([Zs[:,0]-Zs[:,1],Zs[:,1]-Zs[:,2]])
DWRg = np.array([Zg[:,0]-Zg[:,1],Zg[:,1]-Zg[:,2]])
DWRh = np.array([Zh[:,0]-Zh[:,1],Zh[:,1]-Zh[:,2]])

gridsize = (2, 3)
fig = plt.figure(figsize=(12, 6))
ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=2, rowspan=1)
ax2 = plt.subplot2grid(gridsize, (1, 0), colspan=2, rowspan=1)
ax3 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=2)

ax1.plot(avgMassc, DWRc[0,:], lw=2)
ax2.plot(avgMassc, DWRc[1,:], lw=2)
lc = ax3.plot(DWRc[1,:], DWRc[0,:], lw=2)

ax1.plot(avgMassr, DWRr[0,:], lw=2)
ax2.plot(avgMassr, DWRr[1,:], lw=2)
lr = ax3.plot(DWRr[1,:], DWRr[0,:], lw=2)

ax1.plot(avgMassi, DWRi[0,:], lw=2)
ax2.plot(avgMassi, DWRi[1,:], lw=2)
li = ax3.plot(DWRi[1,:], DWRi[0,:], lw=2)

ax1.plot(avgMasss, DWRs[0,:], lw=2)
ax2.plot(avgMasss, DWRs[1,:], lw=2)
ls = ax3.plot(DWRs[1,:], DWRs[0,:], lw=2)

ax1.plot(avgMassg, DWRg[0,:], lw=2)
ax2.plot(avgMassg, DWRg[1,:], lw=2)
lg = ax3.plot(DWRg[1,:], DWRg[0,:], lw=2)

ax1.set_ylim([-2, 15])
ax2.set_ylim([-2, 15])
ax3.set_xlim([-2, 15])
ax3.set_ylim([-2, 15])

ax1.set_xlabel('D$_m$   [mm]')
ax1.set_ylabel('DWR$_{XK_a}$')
ax2.set_xlabel('D$_m$   [mm]')
ax2.set_ylabel('DWR$_{K_aW}$')
ax3.set_xlabel('DWR$_{K_aW}$')
ax3.set_ylabel('DWR$_{XK_a}$')

ax3.legend(handles=[lr[0], li[0], ls[0], lg[0]],
           labels=['rain', 'ice', 'snow', 'graupel'], loc=2)

ax1.grid()
ax2.grid()
ax3.grid()

plt.tight_layout()
plt.savefig('theoretical_DWRs.pdf', dpi=300)
plt.savefig('theoretical_DWRs.png', dpi=300)