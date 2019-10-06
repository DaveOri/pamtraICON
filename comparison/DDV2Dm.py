#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 12:12:41 2019

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
from air_properties import FooteduToit_coeff
from scipy.special import gamma

rho = 1000.0
a = rho*np.pi/6.0
b = 3.0

def mD(D):
  return a*D**b
def Dm(m):
  return (m/a)**(1/b)

av = 114.0137
bv = 0.23437

nu=0.0
mu=1.0/3.0
gam = b*mu
mup = b*nu + b - 1.0

def lam(N, q):
  return (N*gamma((nu+2)/mu)/(q*gamma((nu+1)/mu)))**mu

def A(N, lam):
  return N*mu*lam**((nu+1)/mu)/gamma((nu+1)/mu)

def LAM(lam):
  return lam*a**mu

def N0(A):
  return A*b*a**(nu+1)

def MkN(N0, LAM, k):
  return N0*gamma((mup+k+1)/gam)/(gam*LAM**((mup+k+1)/gam))

def Mkf(A, lam, k):
  return A*gamma((nu+k+1)/mu)/(mu*lam**((nu+k+1)/mu))

def qN2Dm(q, N):
  ll = lam(N, q)
  AA = A(N, ll)
  LL = LAM(ll)
  NN = N0(AA)
  return MkN(NN, LL, 4)/MkN(NN, LL, 3)

def qN2MDVp(q, N):
  ll = lam(N, q)
  AA = A(N, ll)
  return av*Mkf(AA, ll, 2+bv)/Mkf(AA, ll, 2)


def qN2SWp(q, N):
  ll = lam(N, q)
  AA = A(N, ll)
  mdv = qN2MDVp(q, N)
  return av*av*Mkf(AA, ll, 2+2*bv)/Mkf(AA, ll, 2) - mdv

def qN2MDVa(q, N):
  ll = lam(N, q)
  AA = A(N, ll)
  return alpha - beta*Mkf(AA, ll+c, 2)/Mkf(AA, ll, 2)

def qN2SWa(q, N):
  ll = lam(N, q)
  AA = A(N, ll)
  mdv = qN2MDVa(q, N)
  return alpha*alpha + beta*(-2*alpha*Mkf(AA, ll+c, 2) + beta*Mkf(AA, ll+2*c, 2))/Mkf(AA, ll, 2) - mdv

def qN2moments(q, N):
  ll = lam(N, q)
  AA = 1#A(N, ll)
  M2 = Mkf(AA, ll, 2)
  M2_b = Mkf(AA, ll, 2+bv)
  mdv = av*M2_b/M2
  M2b_2 = Mkf(AA, ll, 2*bv+2)
  sw = np.sqrt(av**2*M2b_2/M2 - mdv*mdv)
  M2_3b = Mkf(AA, ll, 2+3*bv)
  sk = (M2_3b*av**3 + M2*mdv**3 -3*M2b_2*mdv*av**2 + 3*M2_b*av*mdv**3)/(M2*sw**3)
  return -mdv, sw, sk

def Matrosov17(DDV): # Ka-W
  if 0.0<=DDV<=1:
    return 0.47+0.49*DDV**0.54
  elif DDV<2.4:
    return 1.338+DDV*(-0.977+DDV*(0.678-DDV*0.079))
  else:
    return np.nan

vMatrosov17 = np.vectorize(Matrosov17)

def Kamil19(DDV): # X-W
  return 0.576+DDV*(0.905+DDV*(-0.779+DDV*(0.451+DDV*(-0.108+DDV*0.009))))

accMins = 5
freq = str(accMins)+'min'
iconfile = 'data/precipitation_icon.h5'
icon = pd.read_hdf(iconfile, key='stat')
icon = icon.reset_index().drop_duplicates(subset='index',
                                          keep='last').set_index('index')
pluviofile = 'data/precipitation_pluvio.h5'
pluvio = pd.read_hdf(pluviofile, key='stat')
mrrfile = 'data/precipitation_mrr.h5'
mrr = pd.read_hdf(mrrfile, key='stat')
mrr = mrr.resample(freq).apply(np.nansum)

ixi = pd.date_range(start='2015-11-11', end='2016-1-4', freq='9s')
icon.reindex(ixi)
icon = icon.resample(freq).apply(np.nansum)

ixp = pd.date_range(start='2015-11-11', end='2016-1-4', freq='1min')
pluvio.reindex(ixp)
pluvio = pluvio.resample(freq).apply(np.nansum)

pamtra = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
                        hydroset='all_hydro', suffix='pamtra_icon.h5',
                        pamtra=True, minhour=6.0,
                        varlist=['Z10', 'Z35', 'Z94', 'W10', 'W35', 'W94',
                                 'T', 'unixtime', 'P', 'RH', 'QNR', 'QR',
                                 'QG', 'QI', 'QS', 'QH', 'QC',
                                 'V10', 'V35', 'V94'])

pamtra['Q'] = pamtra['QR']+pamtra['QG']+pamtra['QI']+pamtra['QS']+pamtra['QH']+pamtra['QC']
pamtra['R/Q'] = pamtra['QR']/pamtra['Q']

#ice = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
#                        hydroset='only_ice', suffix='pamtra_icon.h5', pamtra=True,
#                        varlist=['Z10', 'Z35', 'Z94', 'T',
#                                 'V10', 'V35', 'V94', 'unixtime',
#                                 'W10', 'W35', 'W94'], minhour=6.0)
#
#snow = read_variables(path='/work/develop/pamtraICON/comparison/data/pamtra/',
#                        hydroset='only_snow', suffix='pamtra_icon.h5', pamtra=True,
#                        varlist=['Z10', 'Z35', 'Z94', 'T',
#                                 'V10', 'V35', 'V94', 'unixtime',
#                                 'W10', 'W35', 'W94'], minhour=6.0)

radar = read_variables(path='/work/develop/pamtraICON/comparison/data/radar/',
                       hydroset='', suffix='radar_regrid.h5', minhour=6.0,
                       varlist=['Z10', 'Z35', 'Z94', 'T', 'P', 'RH',
                                'V10avg', 'V35avg', 'V94avg', 'W10', 'W35', 'W94',
                                'unixtime', 'quality_x', 'quality_w'])

radar.unixtime = pd.to_datetime(radar.unixtime.astype(np.int64), unit='s')
pamtra.unixtime = pd.to_datetime(pamtra.unixtime.astype(np.int64), unit='s')

radar['RR'] = (pluvio.resample('1s').nearest().loc[radar.unixtime]*60/accMins).values
pamtra['RR'] = (icon.resample('1s').nearest().loc[pamtra.unixtime]*60/accMins).values

#radarw = slice_data(radar, 'quality_w', maxvalue=8192)
#radarx = slice_data(radar, 'quality_x', maxvalue=8192)

pamtra = slice_data(pamtra, 'Z35', -15.0)
radar = slice_data(radar, 'Z35', -15.0)

logrule = True
density = False
CFAD = True
inverty = True
bins = 100
stats = ['mean', 'median', 'quartile', 'decile']
## High precipitation
minRR = 1.0
maxRR = 91.0
pre = 'HIG'

#mdv, sw, sk = qN2moments(pamtra['QR'], pamtra['QNR'])
#pamtra['MDV'] = mdv
#pamtra['SW'] = sw

f, ((ax11, ax12), (ax21, ax22), (ax31, ax32)) = plt.subplots(3, 2, figsize=(10.5, 9.))
r = hist_and_plot(slice_data(pamtra, 'RR', minvalue=minRR, left=True),
                  'Simulated MDV Ka',
                  yvar='T', xvar='V35',
                  xlabel='MDV   [m/s]', ylabel='T   [deg C]',
                  vminmax=[0.1, 30],
                  xlim=[-10, 0], ylim=[0, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VSD.png',
                  inverty=inverty, figax=(f, ax11), stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radar, 'RR', minvalue=minRR, left=True),
                  'Measured MDV Ka',
                  yvar='T', xvar='V35avg',
                  xlabel='MDV   [m/s]', ylabel='T   [deg C]',
                  vminmax=[0.1, 30],
                  xlim=[-10, 0], ylim=[0, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VSD.png',
                  inverty=inverty, figax=(f, ax12), stats=stats,
                  bins=(r[4], r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radar, 'RR', minvalue=minRR, left=True),
                  'Measured SW Ka',
                  yvar='T', xvar='W35',
                  xlabel='SW   [m/s]',
                  ylabel='T   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 3], ylim=[0, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VSD.png',
                  inverty=inverty, figax=(f, ax22), stats=stats,
                  bins=(np.linspace(0, 3), np.linspace(0, 10)),
                  density=density, CFAD=CFAD)


r = hist_and_plot(slice_data(pamtra, 'RR', minvalue=minRR, left=True),
                  'Simulated SW Ka',
                  yvar='T', xvar='W35',
                  xlabel='SW   [m/s]',
                  ylabel='T   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 3], ylim=[0, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VSD.png',
                  inverty=inverty, figax=(f, ax21), stats=stats,
                  bins=(r[4], r[5]), density=density, CFAD=CFAD)

radar['Dm1'] = vMatrosov17((radar['V94avg']-radar['V35avg'])/FooteduToit_coeff(radar['P'], radar['T']+273.15, radar['RH']))
radar['Dm2'] = Kamil19((radar['V94avg']-radar['V10avg'])/FooteduToit_coeff(radar['P'], radar['T']+273.15, radar['RH']))
pamtra['Dm'] = qN2Dm(pamtra['QR'], pamtra['QNR'])*1000.0
pamtra.loc[pamtra.QNR < 1, 'Dm'] = np.nan
pamtra['Dm1'] = vMatrosov17((pamtra['V94']-pamtra['V35'])/FooteduToit_coeff(pamtra['P'], pamtra['T']+273.15, pamtra['RH']))
pamtra['Dm2'] = Kamil19((pamtra['V94']-pamtra['V10'])/FooteduToit_coeff(pamtra['P'], pamtra['T']+273.15, pamtra['RH']))

radar.loc[radar['Dm1']<0.0]['Dm1'] = np.nan
radar.loc[radar['Dm2']<0.0]['Dm2'] = np.nan


r = hist_and_plot(slice_data(pamtra,
                             'RR', minvalue=minRR, left=True),
                  'Simulated Dm rain',
                  yvar='T', xvar='Dm',
                  xlabel='Dm   [mm]',
                  ylabel='T   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 3], ylim=[0, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VSD.png',
                  inverty=inverty, figax=(f, ax31), stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radar, 'RR', minvalue=minRR, left=True),
                  'Retrieved Dm rain',
                  yvar='T', xvar='Dm2',
                  xlabel='Dm  [mm]',
                  ylabel='T   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 3], ylim=[0, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VSD.png',
                  inverty=inverty, figax=(f, ax32), stats=stats,
                  bins=(r[4], r[5]), density=density, CFAD=CFAD)

f.suptitle('T-MDV  CFADs RR>1 mm/h', fontsize=12, fontweight='heavy', y=0.99)
f.tight_layout(pad=1.5, h_pad=0.5, w_pad=0.5)
f.text(x=0.5, y=0.66, s='T-SW   CFADs', fontsize=12, fontweight='heavy',
       horizontalalignment='center')
f.text(x=0.5, y=0.33, s='T-Dm   CFADs', fontsize=12, fontweight='heavy',
       horizontalalignment='center')
f.savefig(pre+'pamRad_T_VSD.png', dpi=300)

## Low precipitation
minRR = -1.0
maxRR = 1.0
pre = 'LOW'
f, ((ax11, ax12), (ax21, ax22), (ax31, ax32)) = plt.subplots(3, 2, figsize=(10.5, 9.))
r = hist_and_plot(slice_data(pamtra, 'RR', minvalue=minRR, maxvalue=maxRR),
                  'Simulated MDV Ka',
                  yvar='T', xvar='V35',
                  xlabel='MDV   [m/s]', ylabel='T   [deg C]',
                  vminmax=[0.1, 30],
                  xlim=[-10, 0], ylim=[0, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VSD.png',
                  inverty=inverty, figax=(f, ax11), stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radar, 'RR', minvalue=minRR, maxvalue=maxRR),
                  'Measured MDV Ka',
                  yvar='T', xvar='V35avg',
                  xlabel='MDV   [m/s]', ylabel='T   [deg C]',
                  vminmax=[0.1, 30],
                  xlim=[-10, 0], ylim=[0, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VSD.png',
                  inverty=inverty, figax=(f, ax12), stats=stats,
                  bins=(r[4], r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radar, 'RR', minvalue=minRR, maxvalue=maxRR),
                  'Measured SW Ka',
                  yvar='T', xvar='W35',
                  xlabel='SW   [m/s]',
                  ylabel='T   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 3], ylim=[0, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VSD.png',
                  inverty=inverty, figax=(f, ax22), stats=stats,
                  bins=(np.linspace(0, 3), np.linspace(0, 10)),
                  density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(pamtra, 'RR', minvalue=minRR, maxvalue=maxRR),
                  'Simulated SW Ka',
                  yvar='T', xvar='W35',
                  xlabel='SW   [m/s]',
                  ylabel='T   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 3], ylim=[0, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VSD.png',
                  inverty=inverty, figax=(f, ax21), stats=stats,
                  bins=(r[4], r[5]), density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(pamtra,#slice_data(pamtra, 'R/Q', minvalue=0.8),
                  'RR', minvalue=minRR, maxvalue=maxRR),
                  'Simulated Dm rain',
                  yvar='T', xvar='Dm',
                  xlabel='Dm   [mm]',
                  ylabel='T   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 3], ylim=[0, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VSD.png',
                  inverty=inverty, figax=(f, ax31), stats=stats,
                  bins=bins, density=density, CFAD=CFAD)

r = hist_and_plot(slice_data(radar, 'RR', minvalue=minRR, maxvalue=maxRR),
                  'Retrieved Dm rain',
                  yvar='T', xvar='Dm2',
                  xlabel='Dm  [mm]',
                  ylabel='T   [deg C]',
                  vminmax=[0.1, 40],
                  xlim=[0, 3], ylim=[0, 10], lognorm=logrule,
                  savename=pre+'pamRad_T_VSD.png',
                  inverty=inverty, figax=(f, ax32), stats=stats,
                  bins=(r[4], r[5]), density=density, CFAD=CFAD)


f.suptitle('T-MDV  CFADs RR<1 mm/h', fontsize=12, fontweight='heavy', y=0.99)
f.tight_layout(pad=1.5, h_pad=0.5, w_pad=0.5)
f.text(x=0.5, y=0.66, s='T-SW   CFADs', fontsize=12, fontweight='heavy',
       horizontalalignment='center')
f.text(x=0.5, y=0.33, s='T-Dm   CFADs', fontsize=12, fontweight='heavy',
       horizontalalignment='center')
f.savefig(pre+'pamRad_T_VSD.png', dpi=300)