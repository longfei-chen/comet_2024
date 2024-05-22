import copy
import glob
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.modeling import models
from astropy.coordinates import SpectralCoord
from specutils import Spectrum1D
from specutils.manipulation import box_smooth
from specutils.fitting.continuum import fit_continuum
from specutils.fitting import fit_lines
# from uwb_tool import pipeline
from uwb_tool import functions


comet_name = "12P"
obs_date_list = ["20240417", "20240424",
                 "20240510", "20240511", "20240513"]
receiver = "UWB2"
fixed_freq_limit = [1655, 1670]
freq_limit = [1664, 1668]
OH_1667 = 1667.3590 * u.MHz
select_date = obs_date_list[4]


# prepare for the frequency axis
_, freq_array = functions.get_freq_mask_range(receiver, fixed_freq_limit[0], fixed_freq_limit[1])
freq_mask = (freq_array >= freq_limit[0]) & (freq_array <= freq_limit[1])

new_spec_array = copy.deepcopy(freq_array[freq_mask]) * u.MHz
new_spec_coord = SpectralCoord(new_spec_array,
                               doppler_convention="radio",
                               doppler_rest=OH_1667)


# read the doppler corrected ta data files
ta_on_files = "./comet_outputs/{}/{}/product/sourceON/Ta_{}_{}_doppler.npy"
ta_off_files = "./comet_outputs/{}/{}/product/sourceOFF/Ta_{}_{}_doppler.npy"
ta_onoff_files = "./comet_outputs/{}/{}/product/Ta_{}_{}_doppler.npy"

ta_on, ta_off, ta_onoff = {}, {}, {}
for obs_date in obs_date_list:
    ta_on[obs_date] = np.load(ta_on_files.format(comet_name, obs_date, fixed_freq_limit[0], fixed_freq_limit[1]))
    ta_off[obs_date] = np.load(ta_off_files.format(comet_name, obs_date, fixed_freq_limit[0], fixed_freq_limit[1]))
    ta_onoff[obs_date] = np.load(ta_onoff_files.format(comet_name, obs_date, fixed_freq_limit[0], fixed_freq_limit[1]))


ta = ta_onoff[select_date]

sp = Spectrum1D(flux=ta * u.K,
                spectral_axis=new_spec_array,
                velocity_convention="radio")


# Baseline removing
fitted_continuum = fit_continuum(sp)
baseline = fitted_continuum(sp.frequency)
sp_baseline_removed = sp - baseline

# Smoothing
sp_smoothed = box_smooth(sp_baseline_removed, width=5)

# Fitting
init_condition = models.Gaussian1D(-0.4, OH_1667, 0.03*u.MHz)
fit_func = fit_lines(sp_smoothed, init_condition)
sp_fitted = fit_func(sp.frequency)


# plot with respect to frequency
plt.figure(figsize=(10,5))
# plt.plot(sp.frequency, sp.flux, label="raw spectrum")
# plt.plot(sp.frequency, sp_baseline_removed.flux, label="baseline removed")
plt.plot(sp.frequency, sp_smoothed.flux, label="smoothed")
plt.plot(sp.frequency, sp_fitted, label="fitted")


freq_min = sp.frequency[0].value
freq_max = sp.frequency[-1].value
plt.xlim([freq_min, freq_max])
plt.xlabel("Frequency (GHz)", fontsize=16)
plt.ylabel("Ta (K)", fontsize=16)
plt.grid()
plt.title(select_date, fontsize=20)
plt.legend()
plt.show()


# plot with respect to velocity
v_lsr = new_spec_coord.to(u.km / u.s, doppler_convention='radio', doppler_rest=OH_1667)

plt.figure(figsize=(10,5))
plt.plot(v_lsr, sp_smoothed.flux, label="smoothed")
plt.plot(v_lsr, sp_fitted, label="fitted")


v_min = v_lsr[0].value
v_max = v_lsr[-1].value
plt.xlim([v_min, v_max])
plt.xlabel("V$_{lsr}$ (km/s)", fontsize=16)
plt.ylabel("Ta (K)", fontsize=16)
plt.grid()
plt.title(select_date, fontsize=20)
plt.legend()
plt.show()
