#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 11:41:48 2019

@author: dori
"""

import numpy as np

T0 = 273.15
r0 = 1.20386316242422
p0 = 101325.0

def saturation_vapour_pressure(T):
  return 610.78*np.exp(17.27*(T-T0)/(T-T0+237.3))

def RH2q(RH, T, p):
  e = RH*saturation_vapour_pressure(T)/100.0
  return 0.622*e/(p-e)

def dry_air_density(p, T):
  return p/(T*286.9)

def moist_air_density(p, T, RH):
  q = RH2q(RH, T, p)
  return dry_air_density(p, T)*(1+q)/(1+1.609*q)

def FooteduToit_coeff(p, T, RH):
  r = moist_air_density(p, T, RH)
  logr0r = np.log10(r0/r)
  Y = 0.43*logr0r - 0.4*(logr0r)**2.5
  return 10**Y*(1.0+0.0023*(1.1-r/r0)*(T-T0+20.0))