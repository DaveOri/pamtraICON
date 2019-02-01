#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 17:38:02 2019

@author: dori
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from READ import read_prepare
plt.close('all')

hydroset = 'all_hydro'
#hydroset = 'only_snow'
#hydroset = 'only_ice'
hydroset = 'only_graupel_hail'
#hydroset = 'only_liquid'

Nstr = 'QNH'
qstr = 'QH'
hydro = 'Hail'

hydroset='all_hydro'
data = read_prepare(hydroset)

diams = np.linspace(0.00005,0.0229,1000)
mass = lambda d: 0.038*d**2.0
vol = lambda d: (np.pi*d**3.0)/6.0
masses = mass(diams)
volumes = vol(diams)
Dm = lambda lam: ((masses*diams*diams*np.exp(-lam*diams)).sum()) / ((masses*diams*np.exp(-lam*diams)).sum())
D0 = lambda lam: ((volumes*diams*diams*np.exp(-lam*diams)).sum()) / ((volumes*diams*np.exp(-lam*diams)).sum())
vDm = np.vectorize(Dm)
vD0 = np.vectorize(D0)

lam = np.sqrt(0.228*data[Nstr]/data[qstr])
n0 = data[Nstr]*lam**2.0
D0 = vD0(lam)
xlim = [1e-3,5e5]
ylim = [1e-10,1e-2]

# min(data['QI'].values[np.nonzero(data['QI'])])
# min(data['QNI'].values[np.nonzero(data['QNI'])])
# min(data['QS'].values[np.nonzero(data['QS'])])
# min(data['QNS'].values[np.nonzero(data['QNS'])])

plt.close('all')
fig, ax = plt.subplots(1,1)
mapp = ax.scatter(x=data[Nstr], y=data[qstr], s=0.1, c=data['DWRkw'], cmap='jet', vmin=0.0, vmax=15)
plt.colorbar(mapp, ax=ax, label='DWR')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(Nstr + ' [m-3]')
ax.set_ylabel(qstr + ' [kg/kg]')
ax.set_title(hydro + ' moments vs DWR Ka W')
fig.savefig('qNstudy/'+hydro+'_'+hydroset+'MomDWR.png')

fig, ax = plt.subplots(1,1)
mapp = ax.scatter(x=data[Nstr], y=data[qstr], s=0.1, c=D0, cmap='jet', vmin=0.0, vmax=0.01)
plt.colorbar(mapp, ax=ax, label='D0')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(Nstr + ' [m-3]')
ax.set_ylabel(qstr + ' [kg/kg]')
ax.set_title(hydro + ' moments vs D0')
fig.savefig('qNstudy/'+hydro+'_'+hydroset+'MomD0.png')
#plt.close('all')
#fig, ax = plt.subplots(1,1)
#mapp = ax.scatter(x=data['QNS'], y=data['QS'], s=0.1, c=data['DWRkw'], cmap='jet', vmin=0.0, vmax=15)
#plt.colorbar(mapp, ax=ax, label='DWR')
#ax.set_xlim(xlim)
#ax.set_ylim(ylim)
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlabel('QNS [m-3]')
#ax.set_ylabel('QS  [kg/kg]')
#ax.set_title('Snow moments vs DWR Ka W')
#
#fig, ax = plt.subplots(1,1)
#mapp = ax.scatter(x=data['QNR'], y=data['QR'], s=0.1, c=data['DWRkw'], cmap='jet', vmin=0.0, vmax=15)
#plt.colorbar(mapp, ax=ax, label='DWR')
#ax.set_xlim(xlim)
#ax.set_ylim(ylim)
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlabel('QNR [m-3]')
#ax.set_ylabel('QR  [kg/kg]')
#ax.set_title('Rain moments vs DWR Ka W')
#
#fig, ax = plt.subplots(1,1)
#mapp = ax.scatter(x=data['QNC'], y=data['QC'], s=0.1, c=data['DWRkw'], cmap='jet', vmin=0.0, vmax=15)
#plt.colorbar(mapp, ax=ax, label='DWR')
#ax.set_xlim(xlim)
#ax.set_ylim(ylim)
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlabel('QNC [m-3]')
#ax.set_ylabel('QC  [kg/kg]')
#ax.set_title('Cloud droplets moments vs DWR Ka W')
#
#fig, ax = plt.subplots(1,1)
#mapp = ax.scatter(x=data['QNG'], y=data['QG'], s=0.1, c=data['DWRkw'], cmap='jet', vmin=0.0, vmax=15)
#plt.colorbar(mapp, ax=ax, label='DWR')
#ax.set_xlim(xlim)
#ax.set_ylim(ylim)
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlabel('QNG [m-3]')
#ax.set_ylabel('QG  [kg/kg]')
#ax.set_title('Graupel moments vs DWR Ka W')
#
#fig, ax = plt.subplots(1,1)
#mapp = ax.scatter(x=data['QNH'], y=data['QH'], s=0.1, c=data['DWRkw'], cmap='jet', vmin=0.0, vmax=15)
#plt.colorbar(mapp, ax=ax, label='DWR')
#ax.set_xlim(xlim)
#ax.set_ylim(ylim)
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlabel('QNH [m-3]')
#ax.set_ylabel('QH  [kg/kg]')
#ax.set_title('Hail moments vs DWR Ka W')

##plt.close('all')
#fig, ax = plt.subplots(1,1)
#mapp = ax.scatter(x=data['QNS'], y=data['QS'], s=0.1, c=D0s, cmap='jet')
#plt.colorbar(mapp, ax=ax, label='D0')
#ax.set_xlim(xlim)
#ax.set_ylim(ylim)
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlabel('QNS [m-3]')
#ax.set_ylabel('QS  [kg/kg]')
#ax.set_title('Snow moments vs D0')