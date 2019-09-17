#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 17:24:24 2019

@author: dori
"""

import matplotlib.pyplot as plt
import sys
sys.path.append('..')
from READ import slice_data
from READ import read_variables
from statistic import hist_and_plot

pamtra = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='all_hydro', suffix='pamtra_icon.h5',
                        pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T',
                                 'V10', 'V35', 'V94'], minhour=6.0)

ice = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='only_ice', suffix='pamtra_icon.h5',
                        pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T',
                                 'V10', 'V35', 'V94'], minhour=6.0)

snow = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='only_snow', suffix='pamtra_icon.h5',
                        pamtra=True,
                        varlist=['Z10', 'Z35', 'Z94', 'T',
                                 'V10', 'V35', 'V94'], minhour=6.0)

radar = read_variables(path='/work/develop/pamtraICON/comparison/data/radar/',
                       hydroset='', suffix='radar_regrid.h5', minhour=6.0,
                       varlist=['Z10', 'Z35', 'Z94', 'T',
                                'V10avg', 'V35avg', 'V94avg',
                                'quality_x', 'quality_w'])
#rad = slice_data(radar, 'V35avg', minvalue=-2, maxvalue=1)
rad = slice_data(radar, 'Z35', minvalue=-40, maxvalue=30)
#radarw = slice_data(radar, 'quality_w', maxvalue=8192)
#radarx = slice_data(radar, 'quality_x', maxvalue=8192)
#radarxw = slice_data(radarw, 'quality_x', maxvalue=8192)

lognormrule = True
density = True
CFAD = False
bins=100
stats = ['mean', 'median', 'quartile']

pamtra = slice_data(pamtra, 'DWRkw', maxvalue=4)
ice = slice_data(ice, 'DWRkw', maxvalue=4)
snow = slice_data(snow, 'DWRkw', maxvalue=4)
radar = slice_data(radar, 'DWRkw', maxvalue=4)

r = hist_and_plot(slice_data(pamtra, 'T', maxvalue=-2), 'Simulated MDVk - Zk',
                  yvar='V35', xvar='Z35',
              xlabel='Zk   [dBZ]', ylabel='MDV Ka   [m/s]', #vminmax=[0.1, 100],
              xlim=[-40, 30], ylim=[-2, 1], lognorm=lognormrule,
              savename='pamtra_Vk_Zk.png',
              inverty=False, figax=None,
              bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(ice, 'T', maxvalue=-2),
                  'Simulated MDVk - Zk only ICE', yvar='V35', xvar='Z35',
              xlabel='Zk   [dBZ]', ylabel='MDV Ka   [m/s]', #vminmax=[0.1, 100],
              xlim=[-40, 30], ylim=[-2, 1], lognorm=lognormrule,
              savename='pamtra_ICE_Vk_Zk.png',
              inverty=False, figax=None,
              bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(snow, 'T', maxvalue=-2),
                  'Simulated MDVk - Zk only SNOW', yvar='V35', xvar='Z35',
              xlabel='Zk   [dBZ]', ylabel='MDV Ka   [m/s]', #vminmax=[0.1, 100],
              xlim=[-40, 30], ylim=[-2, 1], lognorm=lognormrule,
              savename='pamtra_SNOW_Vk_Zk.png',
              inverty=False, figax=None,
              bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(rad, 'T', maxvalue=-2),
                  'Measured MDVk - Zk', yvar='V35avg', xvar='Z35',
              xlabel='Zk   [dBZ]', ylabel='MDV Ka   [m/s]', #vminmax=[0.1, 100],
              xlim=[-40, 30], ylim=[-2, 1], lognorm=lognormrule,
              savename='radar_Vk_Zk.png',
              inverty=False, figax=None,
              bins=bins, density=density, CFAD=CFAD)

#%% Combined plot for paper

f, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2, figsize=(10.5, 9.))
r = hist_and_plot(slice_data(pamtra, 'T', maxvalue=-2), 'Simulated',
                  yvar='V35', xvar='Z35',
                  xlabel='Zk   [dBZ]', ylabel='MDV Ka   [m/s]', #vminmax=[0.1, 100],
                  xlim=[-40, 30], ylim=[-2, 1], lognorm=lognormrule,
                  savename='pamRad_MDV_Z.png',
                  inverty=False, figax=(f, ax11),
                  bins=bins, density=density, CFAD=CFAD)
ax11.text(0, 0.25, 'All hydrometeors', fontweight='black',
          bbox=dict(facecolor='white', alpha=0.5))

r = hist_and_plot(slice_data(ice, 'T', maxvalue=-2),
                  'Simulated only ICE', yvar='V35', xvar='Z35',
                  xlabel='Zk   [dBZ]', ylabel='MDV Ka   [m/s]', #vminmax=[0.1, 100],
                  xlim=[-40, 30], ylim=[-2, 1], lognorm=lognormrule,
                  savename='pamRad_MDV_Z.png',
                  inverty=False, figax=(f, ax21),
                  bins=bins, density=density, CFAD=CFAD)
ax21.text(-5, 0.25, 'Only ice component', fontweight='black',
          bbox=dict(facecolor='white', alpha=0.5))

r = hist_and_plot(slice_data(snow, 'T', maxvalue=-2),
                  'Simulated only SNOW', yvar='V35', xvar='Z35',
                  xlabel='Zk   [dBZ]', ylabel='MDV Ka   [m/s]', #vminmax=[0.1, 100],
                  xlim=[-40, 30], ylim=[-2, 1], lognorm=lognormrule,
                  savename='pamRad_MDV_Z.png',
                  inverty=False, figax=(f, ax22),
                  bins=bins, density=density, CFAD=CFAD)
ax22.text(-10, 0.25, 'Only snow component', fontweight='black',
          bbox=dict(facecolor='white', alpha=0.5))

r = hist_and_plot(slice_data(rad, 'T', maxvalue=-2), 'Measured',
                  yvar='V35avg', xvar='Z35',
                  xlabel='Zk   [dBZ]', ylabel='MDV Ka   [m/s]', #vminmax=[0.1, 100],
                  xlim=[-40, 30], ylim=[-2, 1], lognorm=lognormrule,
                  savename='pamRad_MDV_Z.png',
                  inverty=False, figax=(f, ax12),
                  bins=bins, density=density, CFAD=CFAD)
ax12.text(0, 0.25, 'Measurements', fontweight='black',
          bbox=dict(facecolor='red', alpha=1.))
autoAxis = ax12.axis()
rec = plt.Rectangle((autoAxis[0], autoAxis[2]),
                    (autoAxis[1]-autoAxis[0]),
                    (autoAxis[3]-autoAxis[2]),
                    fill=False, lw=6, color='r', zorder=0)
rec = ax12.add_patch(rec)
rec.set_clip_on(False)

f.suptitle('MDV - Z', fontsize=12, fontweight='heavy', y=0.99)
f.tight_layout(pad=1.5, h_pad=0.5, w_pad=0.5)
f.savefig('pamRad_MDV_Z.png', dpi=300)