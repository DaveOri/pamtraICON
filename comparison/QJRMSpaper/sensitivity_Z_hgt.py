#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 17:09:37 2019

@author: dori
"""
import sys
sys.path.append('..')
from READ import slice_data
from READ import icon150heights
from READ import read_variables
from statistic import hist_and_plot

pamtra = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='all_hydro', suffix='pamtra_icon.h5',
                        varlist=['Hgt', 'Z10', 'Z35', 'Z94'], minhour=6.0)
radar = read_variables(path='/work/develop/pamtraICON/comparison/data/radar/',
                       hydroset='', suffix='radar.h5',
                       varlist=['Hgt', 'Z10', 'Z35', 'Z94', 'quality_x', 'quality_w'], minhour=6.0)

radarw = slice_data(radar, 'quality_w', maxvalue=8192)
radarx = slice_data(radar, 'quality_x', maxvalue=8192)
#radarxw = slice_data(radarw, 'quality_x', maxvalue=8192)

#%% Z CFAD Height
lognormrule = True
density = True
bins = 50
stats = ['mean', 'median', 'quartile']

r = hist_and_plot(pamtra, 'CFAD Zx hgt', yvar='Hgt', xvar='Z10',
              xlabel='Zx   [dBZ]', ylabel='Hgt   [m]', vminmax=[0.1, 100],
              xlim=[-60, 50], ylim=[0, 12000], lognorm=lognormrule,
              savename='pamtraCFAD_Zx_H.png',
              inverty=False, figax=None, stats=stats,
              bins=(bins, icon150heights[::-1]), density=density, CFAD=True)

r = hist_and_plot(pamtra, 'CFAD Zk hgt', yvar='Hgt', xvar='Z35',
              xlabel='Zk   [dBZ]', ylabel='Hgt   [m]', vminmax=[0.1, 100],
              xlim=[-60, 50], ylim=[0, 12000], lognorm=lognormrule,
              savename='pamtraCFAD_Zk_H.png',
              inverty=False, figax=None, stats=stats,
              bins=(bins, icon150heights[::-1]), density=density, CFAD=True)

r = hist_and_plot(pamtra, 'CFAD Zw hgt', yvar='Hgt', xvar='Z94',
              xlabel='Zw   [dBZ]', ylabel='Hgt   [m]', vminmax=[0.1, 100],
              xlim=[-60, 50], ylim=[0, 12000], lognorm=lognormrule,
              savename='pamtraCFAD_Zw_H.png',
              inverty=False, figax=None, stats=stats,
              bins=(bins, icon150heights[::-1]), density=density, CFAD=True)

r = hist_and_plot(radarx, 'CFAD Zx hgt', yvar='Hgt', xvar='Z10',
              xlabel='Zx   [dBZ]', ylabel='Hgt   [m]', vminmax=[0.1, 100],
              xlim=[-60, 50], ylim=[0, 12000], lognorm=lognormrule,
              savename='radarCFAD_Zx_H.png',
              inverty=False, figax=None, stats=stats,
              bins=(bins, icon150heights[::-1]), density=density, CFAD=True)

r = hist_and_plot(radar, 'CFAD Zk hgt', yvar='Hgt', xvar='Z35',
              xlabel='Zk   [dBZ]', ylabel='Hgt   [m]', vminmax=[0.1, 100],
              xlim=[-60, 50], ylim=[0, 12000], lognorm=lognormrule,
              savename='radarCFAD_Zk_H.png',
              inverty=False, figax=None, stats=stats,
              bins=(bins, icon150heights[::-1]), density=density, CFAD=True)

r = hist_and_plot(radarw, 'CFAD Zw hgt', yvar='Hgt', xvar='Z94',
              xlabel='Zw   [dBZ]', ylabel='Hgt   [m]', vminmax=[0.1, 100],
              xlim=[-60, 50], ylim=[0, 12000], lognorm=lognormrule,
              savename='radarCFAD_Zw_H.png',
              inverty=False, figax=None, stats=stats,
              bins=(bins, icon150heights[::-1]), density=density, CFAD=True)

#%% Combined plot for paper

f, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.5))
#ax1.set_aspect('equal')
#ax2.set_aspect('equal')
r = hist_and_plot(pamtra, 'Simulated', yvar='Hgt', xvar='Z35',
              xlabel='Z$_{K_a}$   [dBZ]', ylabel='Hgt   [m]', vminmax=[0.1, 100],
              xlim=[-60, 50], ylim=[0, 12000], lognorm=lognormrule,
              savename='pamRadCFAD_Zk_H.png',
              inverty=False, figax=(f, ax1), stats=stats+['decile'],
              bins=(bins, icon150heights[::-1]), density=density, CFAD=True)
r = hist_and_plot(radar, 'Measured', yvar='Hgt', xvar='Z35',
              xlabel='Z$_{K_a}$   [dBZ]', ylabel='Hgt   [m]', vminmax=[0.1, 100],
              xlim=[-60, 50], ylim=[0, 12000], lognorm=lognormrule,
              savename='pamRadCFAD_Zk_H.png',
              inverty=False, figax=(f, ax2), stats=stats+['decile'],
              bins=(bins, icon150heights[::-1]), density=density, CFAD=True)
f.suptitle('Reflectivity CFADs', fontsize=12, fontweight='heavy', y=0.99)
f.tight_layout(pad=1.5, h_pad=0.5, w_pad=0.5)
f.savefig('pamRadCFAD_Zk_H.png', dpi=300)

