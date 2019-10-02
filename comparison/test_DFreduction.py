#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 12:49:21 2019

@author: dori
"""

import pandas as pd
import numpy as np

def red(x):
  d = {}
  d['T'] = np.mean(x['T'])
  d['P'] = np.mean(x['P'])
  return pd.Series(d, index=d.keys())

data = {}
data['H'] = [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3]
data['t'] = [0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1]
data['T'] = [1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4]
data['P'] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
DF = pd.DataFrame(data)
g = DF.groupby(['H', 't'])
r = g.apply(red)
