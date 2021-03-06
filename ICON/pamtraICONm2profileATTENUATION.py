from __future__ import division
import pyPamtra
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.dates as md
import time
from datetime import datetime
from sys import argv, path

plt.close('all')

#########################################################################
# PATHS
#########################################################################

# Directory where all the .nc meteograms are stored
ICON_folder = '/data/hdcp/icon/tripex/'
ICON_folder = '/data/inscape/icon/experiments/tripex_220km/'
ICON_folder = '/data/inscape/icon/experiments/tripex_220km/newicon/'

# Folder where your descriptor files are stored (you can use mine, Mario's or the default pamtra)
#descriptor_folder = '/home/mech/workspace/pamtra/descriptorfiles/'
descriptor_folder = '/home/dori/pamtra/descriptorfiles/'

#########################################################################
# FILES
#########################################################################

# Meteogram
#ICON_filename = '1d_vars_DOM02.nc'#_20151124T180000-20151124T200000_c1.nc'
ICON_filename = ICON_folder + 'METEOGRAM_patch001_joyce.nc'

output_nc = 'test_profile_attenuation.nc'
output_Z = 'test_profile_attenuation'
#script, ICON_filename, output_nc, output_Z = argv

# Descriptor file for hydrometeors (Scattering models, m(D), v(D))
#descriptor_filename = 'descriptor_file_2m_liudb.txt'
descriptor_filename = 'descriptor_file_2m_ssrg.txt'

# Directory of pluvios for precipitation comparison
#plufile = '/data/data_hatpro/jue/data/pluvio/201511/pluvio2_jue_20151124.log'
#cols = ['datetime', 'intensity Real Time [mm/h]', 'accumulated RT [mm]',
#        'accumulated Non-RealTime [mm]', 'accumulated total NRT [mm]',
#        'bucket filling level RT [mm]', 'bucket filling level NRT [mm]',
#        'temperature load cell [C]', 'heating status', 'status',
#        'temperature electronics unit [C]', 'supply voltage [V]',
#        'temperature orifice ring rim [C]']
#dtypes = {'datetime':str,'accumulated total NRT [mm]':float}
#def parse_datetime(datestr):
#    if len(datestr) != 18:
#        return np.nan
#    else:
#        t = time.strptime(datestr[:14], '%Y%m%d%H%M%S')
#    return datetime(*t[:6])
#pluvio = pd.read_csv(plufile,sep=';',error_bad_lines=False,index_col=0,#dtype=dtypes,
#                     names=cols,parse_dates=True,date_parser=parse_datetime)
#pluvio = pluvio.reset_index().dropna().set_index('datetime')

plufile = '/data/data_hatpro/jue/data/pluvio/netcdf/1511/pluvio2_jue_20151124.nc'
pluvio = Dataset(plufile, 'r')
pluvars = pluvio.variables
var = 'total_accum_NRT'
#var = 'fill_level_NRT'
TotAccNRT = pluvars[var]
PluvioTime = pluvars['time']
units=PluvioTime.units.split('since')[0]
btPluvio = pd.to_datetime(PluvioTime.units.split('since')[-1])
dtPluvio = pd.to_timedelta([i if i > 0 else np.nan for i in PluvioTime[:]],unit=str(units)[0]) # TODO better than this...
PluvioDateTime = (btPluvio + dtPluvio)#.astype(np.int64)
TotAcc= TotAccNRT[:] - TotAccNRT[0]


#########################################################################
# INIT
#########################################################################
pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptor_folder+descriptor_filename)

# SETTINGS
pam.nmlSet['active'] = True
pam.nmlSet['passive'] = False # Passive is time consuming
pam.nmlSet['radar_attenuation'] = 'disabled'
pam.set['verbose'] = 0 # set verbosity levels
pam.set['pyVerbose'] = 1 # change to 0 if you do not want to see job progress number

#pam.nmlSet["radar_mode"] = "spectrum"

#########################################################################
# READ DATA
ICON_file = Dataset(ICON_filename, mode='r')
#########################################################################
# Akio ICON files give all of the varaibles in a multidimensional array
# ICON.variables["values"] (5-D array)
# ICON files from Tripex are stored in a much more structured way
# LET'S DO IT !!!

vals = ICON_file.variables
Nh = len(vals['height_2'])
Nt = len(vals['time'])
H = np.tile(vals['height_2'][:],(Nt,1))
tt = np.tile(vals['time'][:],(Nh,1)).T
rain = vals['RAIN_GSP'][:] + vals['SNOW_GSP'][:]
print H.shape, tt.shape, rain.shape

# PLOT microphysic quantities (if you want)
#val = vals['QNR'][:]#vals['QNI'][:]+vals['QNS'][:]+vals['QNR'][:]+vals['QNG'][:]+vals['QNH'][:]+vals['QNC'][:] 
#plt.pcolormesh(tt,H,val)
#plt.ylim([0,3000])
#plt.colorbar()
#plt.show()

# First copy the variables out of the netCDF4 file, because I like to waste memory.
# These are 2400x150 matrices of meteograms over Julich
# 2400 are the numer of times
# 150 are the number of levels (height_2)

timeidx = np.arange(4500,4501) 
#timeidx = np.arange(0,Nt)
pamData = dict() # empty dictionary to store pamtra Data

times   = vals['time'] # seconds since 2015-11-24 02:00:03 proplectic gregorian
pamData['hgt'] = np.tile(np.flip(vals['height_2'],0),(len(timeidx),1)) # heights at which fields are defined
units=times.units.split('since')[0]
basetime=pd.to_datetime(times.units.split('since')[-1])
dtimes = pd.to_timedelta(times[timeidx],unit=str(units)[0]) # TODO better than this...
pamData['timestamp'] = (basetime + dtimes).astype(np.int64)

print H.shape, tt.shape, rain.shape, pamData['timestamp'].shape
#plt.figure()
#ax=plt.gca()
#xfmt = md.DateFormatter('%H:%M')
#ax.xaxis.set_major_formatter(xfmt)
#ax.plot(basetime+dtimes,rain,label='icon')
##ax.plot(PluvioDateTime,TotAcc)
#pd.DataFrame(zip(PluvioDateTime,TotAcc),columns=['datetime',var]).dropna().set_index('datetime').plot(ax=ax,label='pluvio')
#ax.legend()
#ax.grid()
#plt.show()

pamData['press']    = np.flip(vals['P'][timeidx],1)    # pressure 
pamData['temp']     = np.flip(vals['T'][timeidx],1)    # temperature
#pamData['wind_u']  = vals['U'][:]    # zonal wind speed
#pamData['wind_v']  = vals['V'][:]    # meridional wind speed
pamData['relhum'] = np.flip(vals['REL_HUM'][timeidx],1)

# Read hydrometeors content
hydro_cmpl = np.zeros((len(timeidx),150,6))
#hydro_cmpl[:,:,#] = vals['QV'][:]   # specific humidity
hydro_cmpl[:,:,0] = np.flip(vals['QC'][timeidx],1)   # specific cloud water content
hydro_cmpl[:,:,1] = np.flip(vals['QI'][timeidx],1)   # specific cloud ice content
hydro_cmpl[:,:,2] = np.flip(vals['QR'][timeidx],1)   # rain mixing ratio
hydro_cmpl[:,:,3] = np.flip(vals['QS'][timeidx],1)   # snow mixing ratio
hydro_cmpl[:,:,4] = np.flip(vals['QG'][timeidx],1)   # graupel mixing ratio
hydro_cmpl[:,:,5] = np.flip(vals['QH'][timeidx],1)   # graupel mixing ratio # TODO report probably error encoding long name ...  should be hail mixing ratio

# Read hydrometeors number concentration
hydro_num_cmpl = np.zeros((len(timeidx),150,6))
hydro_num_cmpl[:,:,0] = np.flip(vals['QNC'][timeidx],1)  # number concentration of cloud water
hydro_num_cmpl[:,:,1] = np.flip(vals['QNI'][timeidx],1)  # number concentration ice
hydro_num_cmpl[:,:,2] = np.flip(vals['QNR'][timeidx],1)  # number concentration droplets
hydro_num_cmpl[:,:,3] = np.flip(vals['QNS'][timeidx],1)  # number concentration snow
hydro_num_cmpl[:,:,4] = np.flip(vals['QNG'][timeidx],1)  # number concentration graupel
hydro_num_cmpl[:,:,5] = np.flip(vals['QNH'][timeidx],1)  # number concentration hail 


pamData["hydro_q"] = hydro_cmpl
pamData["hydro_n"] = hydro_num_cmpl

#pamData["lfrac"] = np.array([1]) ####
#pamData["groundtemp"] = pamData['temp'][149]
#pamData['phalf'] = vals['PHALF'][:]# pressure on the half levels

pam.createProfile(**pamData)

#########################################################################
# RUN
#########################################################################
frequencies = [13.6,35.6,94,220]
cores = 4 # number of parallel cores
pam.runParallelPamtra(np.array(frequencies), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF(output_nc) # SAVE OUTPUT

#########################################################################
# PLOTTING :)
#########################################################################
ZeX = pam.r['Ze'][0,0,:,0,0,0]
ZeK = pam.r['Ze'][0,0,:,1,0,0]
ZeW = pam.r['Ze'][0,0,:,2,0,0]
ZeG = pam.r['Ze'][0,0,:,3,0,0]
ZeX = np.ma.masked_values(ZeX,-9999)
ZeK = np.ma.masked_values(ZeK,-9999)
ZeW = np.ma.masked_values(ZeW,-9999)
ZeG = np.ma.masked_values(ZeG,-9999)

AAX = pam.r['Att_atmo'][0,0,:,0]
AAK = pam.r['Att_atmo'][0,0,:,1]
AAW = pam.r['Att_atmo'][0,0,:,2]
AAG = pam.r['Att_atmo'][0,0,:,3]
AAX = np.ma.masked_values(AAX,-9999)
AAK = np.ma.masked_values(AAK,-9999)
AAW = np.ma.masked_values(AAW,-9999)
AAG = np.ma.masked_values(AAG,-9999)

AHX = pam.r['Att_hydro'][0,0,:,0,0]
AHK = pam.r['Att_hydro'][0,0,:,1,0]
AHW = pam.r['Att_hydro'][0,0,:,2,0]
AHG = pam.r['Att_hydro'][0,0,:,3,0]
AHX = np.ma.masked_values(AHX,-9999)
AHK = np.ma.masked_values(AHK,-9999)
AHW = np.ma.masked_values(AHW,-9999)
AHG = np.ma.masked_values(AHG,-9999)

H = 0.001*pam.p['hgt'][0,0,:]

levels = np.arange(-25,35,1)

figx, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,sharey=True)
figy, ((ay1,ay2),(ay3,ay4)) = plt.subplots(2,2,sharex=True,sharey=True)
figz, ((az1,az2),(az3,az4)) = plt.subplots(2,2,sharex=True,sharey=True)

ax1.plot(ZeX,H,label='dis')
ax2.plot(ZeK,H,label='dis')
ax3.plot(ZeW,H,label='dis')
ax4.plot(ZeG,H,label='dis')

ay1.plot(2.0*AAX.cumsum(),H)
ay2.plot(2.0*AAK.cumsum(),H)
ay3.plot(2.0*AAW.cumsum(),H)
ay4.plot(2.0*AAG.cumsum(),H)

az1.plot(2.0*AHX.cumsum(),H)
az2.plot(2.0*AHK.cumsum(),H)
az3.plot(2.0*AHW.cumsum(),H)
az4.plot(2.0*AHG.cumsum(),H)


pam.nmlSet['radar_attenuation'] = 'top-down'
output_nc = 'test_profile_attenuation_td.nc'
pam.runParallelPamtra(np.array(frequencies), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF(output_nc) # SAVE OUTPUT

ZeX = pam.r['Ze'][0,0,:,0,0,0]
ZeK = pam.r['Ze'][0,0,:,1,0,0]
ZeW = pam.r['Ze'][0,0,:,2,0,0]
ZeG = pam.r['Ze'][0,0,:,3,0,0]
ZeX = np.ma.masked_values(ZeX,-9999)
ZeK = np.ma.masked_values(ZeK,-9999)
ZeW = np.ma.masked_values(ZeW,-9999)
ZeG = np.ma.masked_values(ZeG,-9999)
ax1.plot(ZeX,H,label='t-d')
ax2.plot(ZeK,H,label='t-d')
ax3.plot(ZeW,H,label='t-d')
ax4.plot(ZeG,H,label='t-d')

AAX = pam.r['Att_atmo'][0,0,:,0]
AAK = pam.r['Att_atmo'][0,0,:,1]
AAW = pam.r['Att_atmo'][0,0,:,2]
AAG = pam.r['Att_atmo'][0,0,:,3]
AAX = np.ma.masked_values(AAX,-9999)
AAK = np.ma.masked_values(AAK,-9999)
AAW = np.ma.masked_values(AAW,-9999)
AAG = np.ma.masked_values(AAG,-9999)
AHX = pam.r['Att_hydro'][0,0,:,0,0]
AHK = pam.r['Att_hydro'][0,0,:,1,0]
AHW = pam.r['Att_hydro'][0,0,:,2,0]
AHG = pam.r['Att_hydro'][0,0,:,3,0]
AHX = np.ma.masked_values(AHX,-9999)
AHK = np.ma.masked_values(AHK,-9999)
AHW = np.ma.masked_values(AHW,-9999)
AHG = np.ma.masked_values(AHG,-9999)

ay1.plot(2.0*AAX.cumsum(),H)
ay2.plot(2.0*AAK.cumsum(),H)
ay3.plot(2.0*AAW.cumsum(),H)
ay4.plot(2.0*AAG.cumsum(),H)

az1.plot(2.0*AHX.cumsum(),H)
az2.plot(2.0*AHK.cumsum(),H)
az3.plot(2.0*AHW.cumsum(),H)
az4.plot(2.0*AHG.cumsum(),H)

pam.nmlSet['radar_attenuation'] = 'bottom-up'
output_nc = 'test_profile_attenuation_bu.nc'
pam.runParallelPamtra(np.array(frequencies), pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=cores)
pam.writeResultsToNetCDF(output_nc) # SAVE OUTPUT

ZeX = pam.r['Ze'][0,0,:,0,0,0]
ZeK = pam.r['Ze'][0,0,:,1,0,0]
ZeW = pam.r['Ze'][0,0,:,2,0,0]
ZeG = pam.r['Ze'][0,0,:,3,0,0]
ZeX = np.ma.masked_values(ZeX,-9999)
ZeK = np.ma.masked_values(ZeK,-9999)
ZeW = np.ma.masked_values(ZeW,-9999)
ZeG = np.ma.masked_values(ZeG,-9999)
ax1.plot(ZeX,H,label='b-u')
ax2.plot(ZeK,H,label='b-u')
ax3.plot(ZeW,H,label='b-u')
ax4.plot(ZeG,H,label='b-u')

AAX = pam.r['Att_atmo'][0,0,:,0]
AAK = pam.r['Att_atmo'][0,0,:,1]
AAW = pam.r['Att_atmo'][0,0,:,2]
AAG = pam.r['Att_atmo'][0,0,:,3]
AAX = np.ma.masked_values(AAX,-9999)
AAK = np.ma.masked_values(AAK,-9999)
AAW = np.ma.masked_values(AAW,-9999)
AAG = np.ma.masked_values(AAG,-9999)

AHX = pam.r['Att_hydro'][0,0,:,0,0]
AHK = pam.r['Att_hydro'][0,0,:,1,0]
AHW = pam.r['Att_hydro'][0,0,:,2,0]
AHG = pam.r['Att_hydro'][0,0,:,3,0]
AHX = np.ma.masked_values(AHX,-9999)
AHK = np.ma.masked_values(AHK,-9999)
AHW = np.ma.masked_values(AHW,-9999)
AHG = np.ma.masked_values(AHG,-9999)

ay1.plot(2.0*AAX.cumsum(),H)
ay2.plot(2.0*AAK.cumsum(),H)
ay3.plot(2.0*AAW.cumsum(),H)
ay4.plot(2.0*AAG.cumsum(),H)

az1.plot(2.0*AHX.cumsum(),H)
az2.plot(2.0*AHK.cumsum(),H)
az3.plot(2.0*AHW.cumsum(),H)
az4.plot(2.0*AHG.cumsum(),H)

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

ax1.set_title('X  band')
ax2.set_title('Ka band')
ax3.set_title('W  band')
ax4.set_title('G  band')

ay1.set_title('X  band')
ay2.set_title('Ka band')
ay3.set_title('W  band')
ay4.set_title('G  band')

az1.set_title('X  band')
az2.set_title('Ka band')
az3.set_title('W  band')
az4.set_title('G  band')

ax1.set_ylabel('H    [km]')
ax3.set_ylabel('H    [km]')
ay1.set_ylabel('H    [km]')
ay3.set_ylabel('H    [km]')
az1.set_ylabel('H    [km]')
az3.set_ylabel('H    [km]')

ax3.set_xlabel('Z   [dBZ]')
ax4.set_xlabel('Z   [dBZ]')
ay3.set_xlabel('Atmo Att   [dB]')
ay4.set_xlabel('Atmo Att   [dB]')
az3.set_xlabel('Hydro Att   [dB]')
az4.set_xlabel('Hydro Att   [dB]')

ax1.set_xlim([-30,30])
ax2.set_xlim([-30,30])
ax3.set_xlim([-30,30])
ax4.set_xlim([-30,30])
ax1.set_ylim([0,10])
ax2.set_ylim([0,10])
ax3.set_ylim([0,10])
ax4.set_ylim([0,10])

ay1.set_ylim([0,10])
ay2.set_ylim([0,10])
ay3.set_ylim([0,10])
ay4.set_ylim([0,10])
az1.set_ylim([0,10])
az2.set_ylim([0,10])
az3.set_ylim([0,10])
az4.set_ylim([0,10])



figx.show()
figy.show()
figz.show()

