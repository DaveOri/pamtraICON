#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 18:07:56 2019

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

varlist = ['Z10', 'Z35', 'Z94', 'N10', 'N35', 'N94']

pamtra = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='all_hydro', suffix='pamtra_icon.h5',
                        varlist=varlist, minhour=6.0)

#pamtra[pamtra==9999.0] = np.nan
lognormrule=True
data = 1.0*pamtra#slice_data(pamtra, 'DWRxk', maxvalue=-2)

hxk,xxk,yxk = hist_and_plot(data, 'Nx Nk', yvar='N10', xvar='N35',
                            xlim=[-100, 100], ylim=[-100, 100],
                            xlabel='Nx', ylabel='Nk',
                            lognorm=lognormrule,
                            savename='pamtra3f_Nxk.png', inverty=False, figax=None,
                            bins=100, density=False, CFAD=False)

hkw,xkw,ykw = hist_and_plot(data, 'Nk Nw', yvar='N35', xvar='N94',
                            xlim=[-100, 100], ylim=[-100, 100],
                            xlabel='Nk', ylabel='Nw',
                            lognorm=lognormrule,
                            savename='pamtra3f_Nkw.png', inverty=False, figax=None,
                            bins=100, density=False, CFAD=False)

data = slice_data(pamtra, 'DWRxk', maxvalue=-2)
hxk,xxk,yxk = hist_and_plot(data, 'Nx Nk', yvar='N10', xvar='N35',
                            xlim=[-100, 100], ylim=[-100, 100],
                            xlabel='Nx', ylabel='Nk',
                            lognorm=lognormrule,
                            savename='pamtra3f_Nxk.png', inverty=False, figax=None,
                            bins=100, density=False, CFAD=False)

hkw,xkw,ykw = hist_and_plot(data, 'Nk Nw', yvar='N35', xvar='N94',
                            xlim=[-100, 100], ylim=[-100, 100],
                            xlabel='Nk', ylabel='Nw',
                            lognorm=lognormrule,
                            savename='pamtra3f_Nkw.png', inverty=False, figax=None,
                            bins=100, density=False, CFAD=False)