#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 18:53:47 2019

@author: dori
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('..')
from READ import slice_data
from READ import read_variables
from statistic import hist_and_plot

pamtra = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='all_hydro', suffix='pamtra_icon.h5', pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T'], minhour=6.0)

radar = read_variables(path='/work/develop/pamtraICON/comparison/data/radar/',
                       hydroset='', suffix='radar_regrid.h5', minhour=6.0,
                       varlist=['Z10', 'Z35', 'Z94', 'T',
                                'quality_x', 'quality_w'])

pamtra = slice_data(pamtra, 'Z10', minvalue=-7)
radarw = slice_data(radar, 'quality_w', maxvalue=8192)
radarx = slice_data(radar, 'quality_x', maxvalue=8192)

logrule = True
density = True
CFAD = True
bins = (np.arange(-5, 15, 0.35), np.arange(-45,1,0.6))
vminmax = [0.1, 20]
stats = ['mean', 'median', 'quantile']

r = hist_and_plot(pamtra, 'Simulated CFAD T - DWRxk',
                              yvar='T', xvar='DWRxk',
                              xlabel='DWRxk   [dB]', ylabel='T   [K]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                              savename='pamtra_T_DWRxk.png',
                              inverty=True, figax=None, stats=stats,
                              bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(radarx, 'Measured CFAD T - DWRxk',
                              yvar='T', xvar='DWRxk',
                              xlabel='DWRxk   [dB]', ylabel='T   [K]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                              savename='radar_T_DWRxk.png',
                              inverty=True, figax=None, stats=stats,
                              bins=(r[4],r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(pamtra, 'Simulated CFAD T - DWRkw',
                              yvar='T', xvar='DWRkw',
                              xlabel='DWRkw   [dB]', ylabel='T   [K]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                              savename='pamtra_T_DWRkw.png',
                              inverty=True, figax=None, stats=stats,
                              bins=(r[4],r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(radarw, 'Measured CFAD T - DWRkw',
                              yvar='T', xvar='DWRkw',
                              xlabel='DWRkw   [dB]', ylabel='T   [K]',
                              vminmax=vminmax,
                              xlim=[-5, 15], ylim=[-30, 0], lognorm=logrule,
                              savename='radar_T_DWRkw.png',
                              inverty=True, figax=None, stats=stats,
                              bins=(r[4],r[5]), density=density, CFAD=CFAD)
