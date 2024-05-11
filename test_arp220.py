#!/usr/bin/python3
#-*- encoding: utf-8 -*-
import sys
if sys.version_info.major < 3:
    print("This Code runs under Python 3.")
    exit(0)


import os
import glob
import numpy as np
from astropy.time import Time
import astropy.units as u
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.convolution import convolve, Box1DKernel
import matplotlib.pyplot as plt
from uwb_tool import functions
from uwb_tool import pipeline


obssetup = {}
obssetup["ra"] = 233.738433
obssetup["dec"] = 23.503225
obssetup["mjd"] = 59906.17638888889
obssetup["bandw"] = 1100
obssetup["nchan"] = 1048576

obssetup["sampling_time"] = 1
obssetup["period_on"] = 1
obssetup["period_off"] = 9
obssetup["source_on"] = 300
obssetup["source_switch"] = 30
obssetup["tcal_data"] = "./tcal/20220617-high-cal/0.95-1.75GHz-Th.xlsx"


base_input_path = "./comet_uwb_obs_data"
base_output_path = "./comet_outputs"


target_name = "arp220"

# fits_files = glob.glob("/".join([base_input_path, target_name, "*.fits"]))
# hdu_info = functions.get_fits_info(fits_files[0])
# print(hdu_info)
# hdu_data = functions.get_fits_data(fits_files[0])
# pol_1 = hdu_data[:, :, 0]

rest_frequency = 1667.359   #MHz

freq_mask, freq_array = functions.get_freq_mask_range("UWB2", xmin=1630, xmax=1640)

ta = np.load("/".join([base_output_path, target_name, "product", "sourceON", "Ta_1630_1640_ave.npy"]))
smoothed_ta = convolve(ta, Box1DKernel(12))

baseline = pipeline.medfilter(smoothed_ta, 5, obssetup["bandw"]/obssetup["nchan"])
# baseline = pipeline.polyfit(freq_array, smoothed_ta, order=1)

ta = smoothed_ta - baseline

# methord one
# doppler = pipeline.rest_corr(5377, rest_frequency)
# print(doppler)
# freq_array = pipeline.freq_shift(freq_array, doppler, obssetup["bandw"]/obssetup["nchan"])
# ta = pipeline.ta_shift(ta, doppler, obssetup["bandw"]/obssetup["nchan"])


# methord two
vlsr = pipeline.vlsr_corr(obssetup["ra"], obssetup["dec"], obssetup["mjd"])
print(vlsr)


velo_array = pipeline.freq_to_velo(freq_array, rest_frequency)

#########
plt.plot(velo_array, ta, label="Ta")
# plt.plot(velo_array, baseline, label="smoothed")

plt.xlim([5000,6000])
plt.ylim([-1000,2500])

plt.xlabel("Velo (km/s)")
plt.ylabel("Ta (K)")
plt.grid()
plt.legend()
plt.show()



