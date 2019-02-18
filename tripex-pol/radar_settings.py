# DEFINE LIBRARY OF RADAR PROPRTIES
Joyrad10nml = {'radar_fwhr_beamwidth_deg':1.0,
               'radar_integration_time':2.0,
               'radar_k2':0.93,
               'radar_max_v':78.07291,
               'radar_min_v':-78.07291,
               'radar_nfft':4096,
               'radar_no_ave':10,
               'radar_pnoise0':-48.0,
#               'radar_noise_distance_factor':10.0,
               'frequency':9.6}

Joyrad35nml = {'radar_fwhr_beamwidth_deg':0.6,
               'radar_integration_time':2.0,
               'radar_k2':0.93,
               'radar_max_v':10.56824,
               'radar_min_v':-10.56824,
               'radar_nfft':512,
               'radar_no_ave':38,
               'radar_pnoise0':-64.0,
               'frequency':35.5}

Grarad94nml = {'radar_fwhr_beamwidth_deg':0.5,
               'radar_integration_time':1.0,
               'radar_k2':0.74,
               'radar_max_v':6.8,
               'radar_min_v':-6.8,
               'radar_nfft':512,
               'radar_no_ave':17,
               'radar_pnoise0':-54.0,
               'frequency':94.0}
Default = {}
radarlib = {'Joyrad10':Joyrad10nml,'Joyrad35':Joyrad35nml,'Grarad94':Grarad94nml, 'Default':Default}

# DEFINE DICTIONARY OF HYDROMETOR CONTENT COMBINATIONS
hydrodict = {'all_hydro'        :[1.,1.,1.,1.,1.,1.],
             'no_snow'          :[1.,1.,1.,0.,1.,1.],
             'only_snow'        :[0.,0.,0.,1.,0.,0.],
             'only_ice'         :[0.,1.,0.,0.,0.,0.],
             'only_liquid'      :[1.,0.,1.,0.,0.,0.],
             'only_graupel_hail':[0.,0.,0.,0.,1.,1.]}
