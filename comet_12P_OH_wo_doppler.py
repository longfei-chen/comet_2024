import glob
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.modeling import models
from astropy.coordinates import SpectralCoord
from specutils import Spectrum1D
from specutils.analysis import centroid, gaussian_fwhm, snr_derived
from specutils.manipulation import box_smooth, FluxConservingResampler
from specutils.fitting.continuum import fit_continuum
from specutils.fitting import fit_lines
# from uwb_tool import pipeline
from uwb_tool import functions


# based_data_path = "./"
based_data_path = "F:/comet_2024/"
comet_name = "12P"
# receiver = "UWB2"
# fixed_freq_limit = [1655, 1670]
receiver = "UWB4"
fixed_freq_limit = [3300, 3450]


obs_date_list = ["20240417", "20240424",
                 "20240429", "20240503",
                 "20240510", "20240511", "20240513"]


OH_1665 = 1665.4018 * u.MHz
OH_1667 = 1667.3590 * u.MHz
CH_3264 = 3263.794 * u.MHz
CH_3335 = 3335.481 * u.MHz
CH_3349 = 3349.193 * u.MHz

rest_line = CH_3349
width = 3 # MHz
freq_limit = [np.floor(rest_line.value - width), np.floor(rest_line.value + width)]

select_date = obs_date_list[4]


ta_on_files = f"{based_data_path}" + "comet_outputs/{}/{}/{}/product/sourceON/Ta_{}_{}_ave.npy"
ta_off_files = f"{based_data_path}" + "comet_outputs/{}/{}/{}/product/sourceOFF/Ta_{}_{}_ave.npy"
ta_onoff_files = f"{based_data_path}" + "comet_outputs/{}/{}/{}/product/Ta_{}_{}_ave.npy"

ta_on, ta_off, ta_onoff = {}, {}, {}
for obs_date in obs_date_list:
    ta_on[obs_date] = np.load(ta_on_files.format(comet_name, obs_date, receiver, fixed_freq_limit[0], fixed_freq_limit[1]))
    ta_off[obs_date] = np.load(ta_off_files.format(comet_name, obs_date, receiver, fixed_freq_limit[0], fixed_freq_limit[1]))
    ta_onoff[obs_date] = np.load(ta_onoff_files.format(comet_name, obs_date, receiver, fixed_freq_limit[0], fixed_freq_limit[1]))


_, freq_array = functions.get_freq_mask_range(receiver, fixed_freq_limit[0], fixed_freq_limit[1])
freq_mask = (freq_array >= freq_limit[0]) & (freq_array <= freq_limit[1])

ta = ta_onoff[select_date][freq_mask]

freq_axis = SpectralCoord(freq_array[freq_mask] * u.MHz,
                          doppler_convention="radio",
                          doppler_rest=rest_line)

sp = Spectrum1D(flux=ta * u.K,
                spectral_axis=freq_axis,
                velocity_convention="radio")


# Baseline removing
fitted_continuum = fit_continuum(sp)
baseline = fitted_continuum(sp.frequency)
sp_baseline_removed = sp - baseline

# Smoothing
sp_smoothed = box_smooth(sp_baseline_removed, width=5)

# Fitting
init_condition = models.Gaussian1D(-0.4, rest_line, 0.03*u.MHz)
fit_func = fit_lines(sp_smoothed, init_condition)
sp_fitted = fit_func(sp.frequency)

fit_sp = Spectrum1D(flux=sp_fitted * u.K,
                    spectral_axis=freq_axis)
velo_sp = Spectrum1D(flux=sp_fitted * u.K,
                     spectral_axis=freq_axis.to(u.km/u.s))

rms = functions.get_rms(sp_smoothed.flux)
peak_ta_idx = np.argmax(np.abs(sp_fitted))
peak_ta = sp_fitted[peak_ta_idx]

# print line profile properties
print(f"rms\t\t{rms.to(u.mK):.2f}")
print(f"peak Ta\t\t{peak_ta.to(u.mK):.2f}")
# print(f"SNR: {snr_derived(sp_smoothed, SpectralRegion(1.664*u.GHz, 1.668*u.GHz))}")
print(f"SNR\t\t{np.abs(peak_ta/rms):.2f}")
print(f"FWHM\t\t{gaussian_fwhm(velo_sp):.2f}")
print(f"velocity\t{fit_sp.velocity[peak_ta_idx]:.2f}")


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
v_lsr = freq_axis.to(u.km / u.s, doppler_convention='radio', doppler_rest=rest_line)

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