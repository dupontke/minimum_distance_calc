#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# USAGE:

# PREAMBLE:

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
from plotting_functions import *
from sel_list import *

zeros = np.zeros

min_dist_dat = sys.argv[1]          # 001.040.min_dist_calc.dat
angle_dat = sys.argv[2]             # ../unique_angle_calc/001.040.3evg_gtp_sah_wt_1.unique_angle_calc.ang.dat
system = sys.argv[3]

nSel = len(sel)

bin_size = 0.01
k = 0.001987 # Kcal K^-1 mol^-1
T = 300. # K
kT = k*T
boltz = 2*kT
four_pi = 4*np.pi


# ----------------------------------------
# MAIN PROGRAM:
# ----------------------------------------

# Load in data_file into a numpy array
datalist0 = np.loadtxt(min_dist_dat)
datalist1 = np.loadtxt(angle_dat)

# 2d histogram plot for probability density of angle verse minimum distance data
hist2d(datalist1[:,0],datalist0[:,2],'S143_OG gtp_PA_O3A Angle','S143 gtp_alpha_phosphate Minimum Distance',500,'S143gtp_alpha_OG_PA_O3A_angle.%s' %(system),'joint_probability_density',norm=True,xunits='$\deg$',yunits='$\AA$')
hist2d(datalist1[:,1],datalist0[:,3],'S143_OG gtp_PB_O3B Angle','S143 gtp_beta_phosphate Minimum Distance',500,'S143gtp_beta_OG_PB_O3B_angle.%s' %(system),'joint_probability_density',norm=True,xunits='$\deg$',yunits='$\AA$')

