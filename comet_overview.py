#!/usr/bin/python3
#-*- encoding: utf-8 -*-
import sys
if sys.version_info.major < 3:
    print("This Code runs under Python 3.")
    exit(0)


import os
import glob
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from uwb_tool import functions


obssetup = {}
obssetup["target_ra"] = 0
obssetup["target_dec"] = 0
obssetup["bandw"] = 500
obssetup["nchan"] = 1048576

obssetup["sampling_time"] = 0.1
obssetup["period_on"] = 0.1
obssetup["period_off"] = 9.9
obssetup["source_on"] = 300
obssetup["source_switch"] = 30

obssetup["tcal_data"] = "./tcal/20220617-high-cal/0.95-1.75GHz-Th.xlsx"

tcal_pol_1 = functions.get_tcal(obssetup["tcal_data"])



base_input_path = "./comet_uwb_obs_data"
base_output_path = "./comet_output"

target_name = "arp220"


fits_files = glob.glob("/".join([base_input_path, target_name, "*.fits"]))


hdu_info = functions.get_fits_info(fits_files[0])
print(hdu_info)


freq_mask, freq_array = functions.get_freq_mask_range("UWB-2", xmin=950, xmax=1750)

hdu_data = functions.get_fits_data(fits_files[0])
pol_1 = hdu_data[:, :, 0]


#######
cal = np.average(pol_1[:, :], axis=1)

plt.scatter(np.arange(len(cal)), cal)
plt.xlabel("Time (second)")
plt.ylabel("Power")
plt.grid()
plt.show()


#########
sp = np.average(pol_1, axis=0)
sp = sp[freq_mask]

plt.plot(freq_array, sp)
plt.xlabel("Frequency (MHz)")
plt.ylabel("Power")
plt.grid()
plt.show()



