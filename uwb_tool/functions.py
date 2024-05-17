#!/usr/bin/python3
#-*- encoding: utf-8 -*-

'''
utils functions
'''
import os
import numpy as np
import astropy.io.fits as fits
from scipy.interpolate import interp1d
import pandas as pd

UWB_receiver = {
"UWB1": {"total_range":[0,1100],    "band_width":1100, "recommended_range":[500,950], "cal_file":"0.5-1GHz-Th.xlsx"},
"UWB2": {"total_range":[800,1900],  "band_width":1100, "recommended_range":[950,1750], "cal_file":"0.95-1.75GHz-Th.xlsx"},
"UWB3": {"total_range":[1600,2700], "band_width":1100, "recommended_range":[1750,2550], "cal_file":"1.75-2.55GHz-Th.xlsx"},
"UWB4": {"total_range":[2400,3500], "band_width":1100, "recommended_range":[2550,3300], "cal_file":"2.55-3.3GHz-Th.xlsx"}
}


def prt_info(msg, *args):
    print(">>> " + msg % args)

def get_rms(arr, sigma=None):
    arr_new = np.copy(arr)
    peak = np.nanmax(arr_new)

    rms = np.sqrt(np.nanmean(np.square(arr_new)))
    
    if sigma is None:
        return rms
        
    while peak > sigma*rms:
        arr_new = arr_new[arr_new < sigma*rms]
        rms = np.sqrt(np.nanmean(np.square(arr_new)))
        peak = np.nanmax(arr_new)

    return rms

def shift_arr(arr, offset=0, cval=np.nan):
    res = np.empty_like(arr)

    if offset > 0:
        res[:offset] = cval
        res[offset:] = arr[:-offset]
    elif offset < 0:
        res[offset:] = cval
        res[:offset] = arr[-offset:]
    else:
        res = arr

    return res

def is_dir_exists(path, mkdir=False):
    if os.path.exists(path):
        return True
    elif not mkdir:
        return False

    if path[0] == "/":
        dest = ""
    else:
        dest = "."
    
    for i in path.split("/"):
        if i in ["", "."]:
            continue
        
        dest = "/".join([dest, i])
        
        if not os.path.exists(dest):
            os.mkdir(dest)
    
    return True

def get_tcal(tcal_file, pol_ave=True):
    tcal = pd.read_excel(tcal_file, header=1, usecols="A,C,F")

    tcal = np.array(tcal)

    freq_array = tcal[:, 0]
    
    if pol_ave:
        pol_averaged = np.average(tcal[:, 1:], axis=1)
        tcal_pol = interp1d(freq_array, pol_averaged, kind="linear", bounds_error=False, fill_value=0)
        return tcal_pol

    pol_1 = tcal[:, 1]
    pol_2 = tcal[:, 2]

    tcal_pol_1 = interp1d(freq_array, pol_1, kind="linear", bounds_error=False, fill_value=0)
    tcal_pol_2 = interp1d(freq_array, pol_2, kind="linear", bounds_error=False, fill_value=0)

    return tcal_pol_1, tcal_pol_2


def get_fits_info(fits_file):
    fits_info = {}

    with fits.open(fits_file) as hdul:
        hdu = hdul[1]

        fits_info["nchan"] = hdu.data.field("NCHAN")[0]
        fits_info["chanw"] = hdu.data.field("CHAN_BW")[0]
        fits_info["exposure"] = hdu.data.field("EXPOSURE")[0]
        fits_info["mjd-start"] = hdu.data.field("UTOBS")[0]
        fits_info["mjd-end"] = hdu.data.field("UTOBS")[-1]
        fits_info["utc-start"] = hdu.data.field("DATE-OBS")[0]
        fits_info["utc-end"] = hdu.data.field("DATE-OBS")[-1]

    return fits_info

def get_fits_data(fits_file):
    with fits.open(fits_file) as hdul:
        hdu_data = hdul[1].data.field("DATA")
    
    return hdu_data


def get_freq_mask_range(uwb, xmin=None, xmax=None, nchan=1048576):
    lower, upper = UWB_receiver[uwb]["total_range"]
    
    band_width = UWB_receiver[uwb]["band_width"]
    
    freq_array = np.arange(lower, upper, 1.0*band_width/nchan)

    if xmin is None:
        xmin = lower
    if xmax is None:
        xmax = upper

    freq_mask = (freq_array >= xmin) & (freq_array <= xmax)

    return freq_mask, freq_array[freq_mask]


def clean_restfreqlist(restfreqs, sel_freq_range_list):
    new_restfreq_list = []

    prev_line = 0
    for line in restfreqs:
        #two transitions should not be close to 0.7 MHz
        if line - prev_line < 0.7:
            continue
        
        prev_line = line

        #is this transition line in the selected frequency range list
        for freq_range in sel_freq_range_list:
            fmin, fmax = freq_range

            if line < fmin + 0.7 or line > fmax-0.7:
                continue
            
            new_restfreq_list.append(line)

    return new_restfreq_list
