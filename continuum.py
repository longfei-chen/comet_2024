import os
import sys
import copy
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.modeling import models
from astropy.coordinates import SpectralCoord
from specutils import Spectrum1D, SpectralRegion
from specutils.analysis import centroid, gaussian_fwhm, snr_derived
from specutils.manipulation import box_smooth
from specutils.fitting.continuum import fit_continuum
from specutils.fitting import fit_lines

from uwb_tool import pipeline
from uwb_tool import functions
from uwb_tool import plotting
from uwb_tool import linedb
from uwb_tool import molspec
import utils


load_path = "F:/comet_2024/"
comet_name = "12P"

obs_date_list = ["20240417", "20240424", 
                 "20240429", "20240503", 
                 "20240510", "20240511", "20240513"]

ta_on_file = f"{load_path}" + "comet_outputs/{}/{}/{}/product/sourceON/Ta_{}_{}_doppler.npy"
ta_off_file = f"{load_path}" + "comet_outputs/{}/{}/{}/product/sourceOFF/Ta_{}_{}_doppler.npy"
ta_onoff_file = f"{load_path}" + "comet_outputs/{}/{}/{}/product/Ta_{}_{}_doppler.npy"


full_band_ta = []
full_band_freq = []

for eachobs in obs_date_list:
    if eachobs not in ["20240429"]:
        continue

    for receiver in utils.band_freq_range_dict.keys():
        frequency_intervals = len(utils.band_freq_range_dict[receiver])
        for i in range(0, frequency_intervals-1):
            freq_lower, freq_upper = utils.band_freq_range_dict[receiver][i:i+2]

            ta = np.load(ta_onoff_file.format(comet_name, eachobs, receiver, freq_lower, freq_upper))

            _, freq_array = functions.get_freq_mask_range(receiver, freq_lower, freq_upper)

            full_band_ta = np.concatenate((full_band_ta, ta))
            full_band_freq = np.concatenate((full_band_freq, freq_array))


spec_coord = SpectralCoord(full_band_freq * u.MHz)
sp = Spectrum1D(flux=full_band_ta * u.K, spectral_axis=spec_coord)

plt.plot(sp.frequency, sp.flux)

plt.show()


