#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 17:38:02 2019

@author: dori
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.close('all')

data = pd.read_hdf('data_pamtra_icon.h5',key='stat')
data[data==-9999] = np.nan
dataZ = data.dropna()

dataZ['DWRxk'] = dataZ.loc[:,'Z10'] - dataZ.loc[:,'Z35']
dataZ['DWRkw'] = dataZ.loc[:,'Z35'] - dataZ.loc[:,'Z94']

diams = np.linspace(0.00005,0.0229,1000)
mass = lambda d: 0.038*d**2.0
vol = lambda d: (np.pi*d**3.0)/6.0
masses = mass(diams)
volumes = vol(diams)
Dm = lambda lam: ((masses*diams*diams*np.exp(-lam*diams)).sum()) / ((masses*diams*np.exp(-lam*diams)).sum())
D0 = lambda lam: ((volumes*diams*diams*np.exp(-lam*diams)).sum()) / ((volumes*diams*np.exp(-lam*diams)).sum())
vDm = np.vectorize(Dm)
vD0 = np.vectorize(D0)

lams = np.sqrt(0.228*dataZ['QNS']/dataZ['QS'])
n0s = dataZ['QNS']*lams**2.0
D0s = vD0(lams)
#Dms = vDm(lams)

lami = np.sqrt(0.228*dataZ['QNI']/dataZ['QI'])
n0i = dataZ['QNI']*lami**2.0
D0i = vD0(lami)
xlim = [1e-30,5e5]
ylim = [1e-30,1e-2]

min(dataZ['QI'].values[np.nonzero(dataZ['QI'])])
min(dataZ['QNI'].values[np.nonzero(dataZ['QNI'])])
min(dataZ['QS'].values[np.nonzero(dataZ['QS'])])
min(dataZ['QNS'].values[np.nonzero(dataZ['QNS'])])

plt.close('all')
fig, ax = plt.subplots(1,1)
mapp = ax.scatter(x=dataZ['QNI'], y=dataZ['QI'], s=0.1, c=dataZ['DWRkw'], cmap='jet', vmin=0.0, vmax=15)
plt.colorbar(mapp, ax=ax, label='DWR')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('QNI [m-3]')
ax.set_ylabel('Qi  [kg/kg]')
ax.set_title('Ice moments vs DWR Ka W')

#plt.close('all')
fig, ax = plt.subplots(1,1)
mapp = ax.scatter(x=dataZ['QNS'], y=dataZ['QS'], s=0.1, c=dataZ['DWRkw'], cmap='jet', vmin=0.0, vmax=15)
plt.colorbar(mapp, ax=ax, label='DWR')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('QNS [m-3]')
ax.set_ylabel('QS  [kg/kg]')
ax.set_title('Snow moments vs DWR Ka W')

fig, ax = plt.subplots(1,1)
mapp = ax.scatter(x=dataZ['QNR'], y=dataZ['QR'], s=0.1, c=dataZ['DWRkw'], cmap='jet', vmin=0.0, vmax=15)
plt.colorbar(mapp, ax=ax, label='DWR')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('QNR [m-3]')
ax.set_ylabel('QR  [kg/kg]')
ax.set_title('Rain moments vs DWR Ka W')

fig, ax = plt.subplots(1,1)
mapp = ax.scatter(x=dataZ['QNC'], y=dataZ['QC'], s=0.1, c=dataZ['DWRkw'], cmap='jet', vmin=0.0, vmax=15)
plt.colorbar(mapp, ax=ax, label='DWR')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('QNC [m-3]')
ax.set_ylabel('QC  [kg/kg]')
ax.set_title('Cloud droplets moments vs DWR Ka W')

fig, ax = plt.subplots(1,1)
mapp = ax.scatter(x=dataZ['QNG'], y=dataZ['QG'], s=0.1, c=dataZ['DWRkw'], cmap='jet', vmin=0.0, vmax=15)
plt.colorbar(mapp, ax=ax, label='DWR')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('QNG [m-3]')
ax.set_ylabel('QG  [kg/kg]')
ax.set_title('Graupel moments vs DWR Ka W')

fig, ax = plt.subplots(1,1)
mapp = ax.scatter(x=dataZ['QNH'], y=dataZ['QH'], s=0.1, c=dataZ['DWRkw'], cmap='jet', vmin=0.0, vmax=15)
plt.colorbar(mapp, ax=ax, label='DWR')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('QNH [m-3]')
ax.set_ylabel('QH  [kg/kg]')
ax.set_title('Hail moments vs DWR Ka W')

plt.close('all')
fig, ax = plt.subplots(1,1)
mapp = ax.scatter(x=dataZ['QNI'], y=dataZ['QI'], s=0.1, c=D0i, cmap='jet', vmin=0.0, vmax=0.01)
plt.colorbar(mapp, ax=ax, label='D0')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('QNI [m-3]')
ax.set_ylabel('QI  [kg/kg]')
ax.set_title('Ice moments vs D0')

#plt.close('all')
fig, ax = plt.subplots(1,1)
mapp = ax.scatter(x=dataZ['QNS'], y=dataZ['QS'], s=0.1, c=D0s, cmap='jet')
plt.colorbar(mapp, ax=ax, label='D0')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('QNS [m-3]')
ax.set_ylabel('QS  [kg/kg]')
ax.set_title('Snow moments vs D0')