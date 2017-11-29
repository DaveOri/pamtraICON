# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 17:34:01 2017

@author: dori
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

icon_folder = '/data/inscape/icon/experiments/'

mid_lat_filename   = icon_folder + 'tripex_220km/newicon/METEOGRAM_patch001_joyce.nc'
mix_phase_filename = icon_folder + 'acloud/RF06-2705/postproc/rf06_sat_all.nc'

mid_lat_data = Dataset(mid_lat_filename)
mix_phase_data = Dataset(mix_phase_filename)

ml_dt_dim = mid_lat_data.dimensions
ml_dt_var = mid_lat_data.variables

mp_dt_dim = mix_phase_data.dimensions
mp_dt_var = mix_phase_data.variables

def plot_variable(variables, dimensions, variable_name):
    var = variables[variable_name]
    var_dim_labels = var.dimensions
    