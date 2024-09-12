import os
import sys
import copy
import numpy as np
from scipy.optimize import curve_fit
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


# load_path = "./"
load_path = "F:/comet_2024/"
comet_name = "12P"

obs_date_list = ["20240417", "20240424", 
                 "20240429", "20240503", 
                 "20240510", "20240511", "20240513"]

def baseline_func(x, a, b):
    return a * x + b

def power_law_func(x, a, b):
    return a * np.power(x, b)


def get_full_band_sp(obs_date, freq_range=None):
    # Ta to Tmb
    # mb_efficiency = 0.85 # Ta/Tmb

    # flux density (Jy) to Ta for FAST
    if freq_range is not None:
        gain = utils.get_gain(np.mean(freq_range))
    else:
        gain = 1.0

    full_band_ta = []
    full_band_freq = []

    ta_on_file = f"{load_path}" + "comet_outputs/{}/{}/{}/product/sourceON/Ta_{}_{}_ave.npy"
    ta_off_file = f"{load_path}" + "comet_outputs/{}/{}/{}/product/sourceOFF/Ta_{}_{}_ave.npy"
    ta_onoff_file = f"{load_path}" + "comet_outputs/{}/{}/{}/product/Ta_{}_{}_ave.npy"

    for receiver in utils.band_freq_range_dict.keys():
        frequency_intervals = len(utils.band_freq_range_dict[receiver])
        for i in range(0, frequency_intervals-1):
            freq_lower, freq_upper = utils.band_freq_range_dict[receiver][i:i+2]

            ta = np.load(ta_on_file.format(comet_name, obs_date, receiver, freq_lower, freq_upper))
            # ta = ta.mean(axis=0)

            _, freq_array = functions.get_freq_mask_range(receiver, freq_lower, freq_upper)

            full_band_ta = np.concatenate((full_band_ta, ta))
            full_band_freq = np.concatenate((full_band_freq, freq_array))
    
    if freq_range is not None:
        freq_range_index = np.where((full_band_freq >= freq_range[0]) & (full_band_freq <= freq_range[1]))

        full_band_freq = full_band_freq[freq_range_index]
        full_band_ta = full_band_ta[freq_range_index]

    spec_coord = SpectralCoord(full_band_freq * u.MHz)
    sp = Spectrum1D(flux=full_band_ta/gain * u.K, spectral_axis=spec_coord)

    return sp


rfi_free_bands = [
[2300,2350], [2375,2399], [2525,2550], [2600,2650], [2700,2760], 
[2775,2800], [2850,2925], [3000,3040], [3150,3300]]

idx = 0
sp_20240417 = get_full_band_sp("20240417", rfi_free_bands[idx])
sp_20240429 = get_full_band_sp("20240429", rfi_free_bands[idx])

sp_20240424 = get_full_band_sp("20240424", rfi_free_bands[idx])
sp_20240503 = get_full_band_sp("20240503", rfi_free_bands[idx])

param_0424, _ = curve_fit(baseline_func, sp_20240424.frequency, sp_20240424.flux)
param_0503, _ = curve_fit(baseline_func, sp_20240503.frequency, sp_20240503.flux)

continuum_0424 = baseline_func(sp_20240424.frequency.value, *param_0424)
continuum_0503 = baseline_func(sp_20240503.frequency.value, *param_0503)

continuum = continuum_0424 - continuum_0503
opt_cont, _ = curve_fit(power_law_func, sp_20240424.frequency, continuum, maxfev=10000)


# sp_20240510 = get_full_band_sp("20240510", rfi_free_bands[idx])
# sp_20240511 = get_full_band_sp("20240511", rfi_free_bands[idx])
# sp_20240513 = get_full_band_sp("20240513", rfi_free_bands[idx])



fig, axs = plt.subplots(1, 2, figsize=(10,4), sharex=True, sharey=True)


# sp1 = sp_20240417 - sp_20240429
# sp2 = sp_20240424 - sp_20240503
# fig.suptitle("(comet + sky_background) - sky_background")
# axs[0].plot(sp1.frequency, sp1.flux, label="20240417 - 20240429")
# axs[1].plot(sp2.frequency, sp2.flux, label="20240424 - 20240503")

fig.suptitle("red: comet + sky_background\nblue: sky_background")
axs[0].plot(sp_20240417.frequency, sp_20240417.flux, label="20240417", color="red")
axs[0].plot(sp_20240429.frequency, sp_20240429.flux, label="20240429", color="blue")

axs[1].plot(sp_20240424.frequency, sp_20240424.flux, label="20240424", color="red")
axs[1].plot(sp_20240503.frequency, sp_20240503.flux, label="20240503", color="blue")
axs[1].plot(sp_20240424.frequency, continuum_0424, color="red")
axs[1].plot(sp_20240503.frequency, continuum_0503, color="blue")


axs[0].set_xlabel("Frequency (GHz)")
axs[0].set_ylabel("Flux (Jy)")
axs[0].grid()
axs[0].legend()

axs[1].set_xlabel("Frequency (GHz)")
axs[1].set_ylabel("Flux (Jy)")
axs[1].grid()
axs[1].legend()
plt.show()


# continuum-subtracted spectrum
rms = functions.get_rms(sp_20240424.flux.value[0:18000]-continuum_0424[0:18000], sigma=3)
print(f"rms: {rms:.2f} Jy")
print(f"alpha: {opt_cont[1]}")
# plt.plot(sp_20240424.frequency, sp_20240424.flux.value-continuum_0424, label="20240424")
# plt.plot(sp_20240503.frequency, sp_20240503.flux.value-continuum_0503, label="20240503")

plt.plot(sp_20240424.frequency, continuum, label="continuum")
plt.plot(sp_20240424.frequency, power_law_func(sp_20240424.frequency.value, *opt_cont))

plt.xlabel("Frequency (GHz)")
plt.ylabel("Flux (Jy)")
plt.grid()
plt.legend()
plt.show()


# plt.plot(sp_20240510.frequency, sp_20240510.flux, label="20240510")
# plt.plot(sp_20240511.frequency, sp_20240511.flux, label="20240511")
# plt.plot(sp_20240513.frequency, sp_20240513.flux, label="20240513")

# plt.xlabel("Frequency (GHz)")
# plt.ylabel("Flux (Jy)")
# plt.grid()
# plt.legend()
# plt.show()
