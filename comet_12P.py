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
obssetup["receiver"] = "UWB2"
obssetup["bandw"] = functions.UWB_receiver[obssetup["receiver"]]["band_width"]
obssetup["nchan"] = 1048576

obssetup["sampling_time"] = 1
obssetup["period_on"] = 1
obssetup["period_off"] = 9
obssetup["source_on"] = 300
obssetup["source_switch"] = 30
obssetup["tcal_data"] = "./tcal/20220617-high-cal/0.95-1.75GHz-Th.xlsx"

comet_name = "12P"
obsdate = "20240510"

base_output_path = f"./comet_outputs/{comet_name}/{obsdate}"


xmin, xmax = 1418, 1422
freq_mask, freq_array = functions.get_freq_mask_range(obssetup["receiver"], xmin=xmin, xmax=xmax)

HI_on = np.load(f"{base_output_path}/product/sourceON/Ta_{xmin}_{xmax}_ave.npy")
HI_off = np.load(f"{base_output_path}/product/sourceOFF/Ta_{xmin}_{xmax}_ave.npy")
HI = np.load(f"{base_output_path}/product/Ta_{xmin}_{xmax}_ave.npy")


plt.plot(freq_array, HI)

plt.xlabel("Freqency (MHz)")
plt.ylabel("Ta (K)")
plt.grid()
plt.show()


xmin, xmax = 1656, 1668
freq_mask, freq_array = functions.get_freq_mask_range(obssetup["receiver"], xmin=xmin, xmax=xmax)

OH_on = np.load(f"{base_output_path}/product/sourceON/Ta_{xmin}_{xmax}_ave.npy")
OH_off = np.load(f"{base_output_path}/product/sourceOFF/Ta_{xmin}_{xmax}_ave.npy")
OH = np.load(f"{base_output_path}/product/Ta_{xmin}_{xmax}_ave.npy")

plt.plot(freq_array, OH)

plt.xlabel("Freqency (MHz)")
plt.ylabel("Ta (K)")
plt.grid()
plt.show()


