#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 17:33:42 2019

@author: dori
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import sys
sys.path.append('..')
from READ import slice_data
from READ import read_variables
from statistic import hist_and_plot

varlist = ['Z10', 'Z35', 'Z94', 'QI', 'QS', 'QNI', 'QNS']

pamtra = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='all_hydro', suffix='pamtra_icon.h5',
                        varlist=varlist, minhour=6.0)

#pamtra[pamtra==9999.0] = np.nan
lognormrule=True
pamtra['QNI'] = np.log10(pamtra['QNI'])
pamtra['QNS'] = np.log10(pamtra['QNS'])
pamtra['QI'] = np.log10(pamtra['QI'])
pamtra['QS'] = np.log10(pamtra['QS'])
pamtra[~np.isfinite(pamtra)] = np.nan

data = 1.0*pamtra#slice_data(pamtra, 'DWRxk', maxvalue=-2)
#data = slice_data(pamtra, 'DWRxk', maxvalue=-2)

h,x,y = hist_and_plot(data, '3f plot all', yvar='QI', xvar='QNI',
                      xlim=[-300, 0], ylim=[-300, 0],
                      xlabel='QI', ylabel='QNI',
                      lognorm=lognormrule,
                      savename='pamtra3f_ice.png', inverty=False, figax=None,
                      bins=100, density=False, CFAD=False)

h,x,y = hist_and_plot(data, '3f plot all', yvar='QS', xvar='QNS',
                      xlim=[-300, 0], ylim=[-300, 0],
                      xlabel='QS', ylabel='QNS',
                      lognorm=lognormrule,
                      savename='pamtra3f_snow.png', inverty=False, figax=None,
                      bins=100, density=False, CFAD=False)

data = slice_data(pamtra, 'DWRxk', maxvalue=-2)

h,x,y = hist_and_plot(data, '3f plot all', yvar='QI', xvar='QNI',
                      xlim=[-300, 0], ylim=[-300, 0],
                      xlabel='QI', ylabel='QNI',
                      lognorm=lognormrule,
                      savename='pamtra3f_ice.png', inverty=False, figax=None,
                      bins=100, density=False, CFAD=False)

h,x,y = hist_and_plot(data, '3f plot all', yvar='QS', xvar='QNS',
                      xlim=[-300, 0], ylim=[-300, 0],
                      xlabel='QS', ylabel='QNS',
                      lognorm=lognormrule,
                      savename='pamtra3f_snow.png', inverty=False, figax=None,
                      bins=100, density=False, CFAD=False)