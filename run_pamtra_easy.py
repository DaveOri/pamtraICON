import sys
import pyPamtra # main pamtra module
import os # import operating system module

# os.environ['PAMTRA_DATADIR'] = '/net/sever/mech/pamtra/data/'

pam = pyPamtra.pyPamtra() # basic empty pyPamtra object with default settings
pam.df.addHydrometeor(("ice", -99., -1, 917., 130., 3.0, 0.684, 2., 3, 1, "mono_cosmo_ice", -99., -99., -99., -99., -99., -99., "mie-sphere", "heymsfield10_particles",0.0))

pamData = dict()
pamData["temp"] = [275.,274.,273.]
pamData["relhum"] = [10.,90.,90.]
pamData["hgt"] = [400.,1500.,2500.]
pamData["press"] = [90000.,80000.,70000.]

pam.createProfile(**pamData) # create a pamtra profile.

pam.runPamtra(89.0)
