#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 16:35:02 2019

@author: dori
"""

import sys
sys.path.append('..')
from READ import slice_data
from READ import read_variables
from statistic import hist_and_plot
import matplotlib.pyplot as plt
import netCDF4

pamtra = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='all_hydro', suffix='pamtra_icon.h5', pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T',
                                 'W10', 'W35', 'W94'], minhour=6.0)

ice = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='only_ice', suffix='pamtra_icon.h5', pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T',
                                 'W10', 'W35', 'W94'], minhour=6.0)

snow = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='only_snow', suffix='pamtra_icon.h5', pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T',
                                 'W10', 'W35', 'W94'], minhour=6.0)

radar = read_variables(path='/work/develop/pamtraICON/comparison/data/radar/',
                       hydroset='', suffix='radar_regrid.h5', minhour=6.0,
                       varlist=['Z10', 'Z35', 'Z94', 'T',
                                'W10', 'W35', 'W94', 'quality_x', 'quality_w'])

radarw = slice_data(radar, 'quality_w', maxvalue=8192)
radarx = slice_data(radar, 'quality_x', maxvalue=8192)

logrule = True
density = False
CFAD = True
inverty = True
bins = 100
stats = ['mean', 'median', 'quartile']

r = hist_and_plot(pamtra, 'Simulated CFAD   T - SW',
                              yvar='T', xvar='W35',
                              xlabel='Spectral Width Ka   [m/s]', ylabel='T   [K]',
                              vminmax=[0.1, 40],
                              xlim=[0, 3], ylim=[-30, 10], lognorm=logrule,
                              savename='pamtra_T_SWk.png',
                              inverty=inverty, figax=None, stats=stats,
                              bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(radar, 'Measured CFAD   T- SW',
                              yvar='T', xvar='W35',
                              xlabel='Spectral Width Ka   [m/s]', ylabel='T   [K]',
                              vminmax=[0.1, 40],
                              xlim=[0, 3], ylim=[-30, 10], lognorm=logrule,
                              savename='radar_T_SWk.png',
                              inverty=inverty, figax=None, stats=stats,
                              bins=bins, density=density, CFAD=CFAD)