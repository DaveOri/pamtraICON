#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 19:48:03 2019

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
                                 'V10', 'V35', 'V94'], minhour=6.0)

ice = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='only_ice', suffix='pamtra_icon.h5', pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T',
                                 'V10', 'V35', 'V94'], minhour=6.0)

snow = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='only_snow', suffix='pamtra_icon.h5', pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T',
                                 'V10', 'V35', 'V94'], minhour=6.0)

radar = read_variables(path='/work/develop/pamtraICON/comparison/data/radar/',
                       hydroset='', suffix='radar.h5', minhour=6.0,
                       varlist=['Z10', 'Z35', 'Z94', 'T',
                                'V10avg', 'V35avg', 'V94avg', 'quality_x', 'quality_w'])

radarw = slice_data(radar, 'quality_w', maxvalue=8192)
radarx = slice_data(radar, 'quality_x', maxvalue=8192)

logrule = True
density = False
CFAD = True
inverty = True
bins = 100

hpw, xpw, ypw = hist_and_plot(pamtra, 'Simulated CFAD   T - MDV',
                              yvar='T', xvar='V35',
                              xlabel='MDV ka   [m/s]', ylabel='T   [K]',
                              vminmax=[0.1, 30],
                              xlim=[-5, 1], ylim=[-30, 10], lognorm=logrule,
                              savename='pamtra_T_Vk.png',
                              inverty=inverty, figax=None,
                              bins=bins, density=density, CFAD=CFAD)

hrx, xrx, yrx = hist_and_plot(radar, 'Measured CFAD   T- MDV',
                              yvar='T', xvar='V35avg',
                              xlabel='MDV ka   [m/s]', ylabel='T   [K]',
                              vminmax=[0.1, 30],
                              xlim=[-5, 1], ylim=[-30, 10], lognorm=logrule,
                              savename='radar_T_Vk.png',
                              inverty=inverty, figax=None,
                              bins=bins, density=density, CFAD=CFAD)