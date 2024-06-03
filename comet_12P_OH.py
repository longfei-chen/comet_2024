import copy
# import glob
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


comet_name = "12P"
obs_date_list = ["20240417", "20240424",
                 "20240510", "20240511", "20240513"]
receiver = "UWB2"
fixed_freq_limit = [1655, 1670]

OH_1667 = 1667.3590 * u.MHz
freq_limit = [1664, 1668]

OH_1665 = 1665.4018 * u.MHz
# freq_limit = [1664, 1667]

select_date = obs_date_list[0]


# prepare for the frequency axis
_, freq_array = functions.get_freq_mask_range(receiver, fixed_freq_limit[0], fixed_freq_limit[1])
freq_mask = (freq_array >= freq_limit[0]) & (freq_array <= freq_limit[1])

new_spec_array = copy.deepcopy(freq_array[freq_mask]) * u.MHz
new_spec_coord = SpectralCoord(new_spec_array,
                               doppler_convention="radio",
                               doppler_rest=OH_1665)


# read the doppler corrected ta data files
ta_on_files = "./comet_outputs/{}/{}/product/sourceON/Ta_{}_{}_doppler.npy"
ta_off_files = "./comet_outputs/{}/{}/product/sourceOFF/Ta_{}_{}_doppler.npy"
ta_onoff_files = "./comet_outputs/{}/{}/product/Ta_{}_{}_doppler.npy"

ta_on, ta_off, ta_onoff = {}, {}, {}
for obs_date in obs_date_list:
    ta_on[obs_date] = np.load(ta_on_files.format(comet_name, obs_date, fixed_freq_limit[0], fixed_freq_limit[1]))
    ta_off[obs_date] = np.load(ta_off_files.format(comet_name, obs_date, fixed_freq_limit[0], fixed_freq_limit[1]))
    ta_onoff[obs_date] = np.load(ta_onoff_files.format(comet_name, obs_date, fixed_freq_limit[0], fixed_freq_limit[1]))


if select_date in ["20240511", "20240513"]:
    ta = (ta_onoff["20240510"] + ta_onoff["20240511"] + ta_onoff["20240513"])/3
else:
    ta = ta_onoff[select_date]

sp = Spectrum1D(flux=ta * u.K,
                spectral_axis=new_spec_coord,
                velocity_convention="radio")


# Baseline removing
fitted_continuum = fit_continuum(sp)
baseline = fitted_continuum(sp.frequency)
sp_baseline_removed = sp - baseline


# Smoothing
sp_smoothed = box_smooth(sp_baseline_removed, width=5)

# standing wave removing
if select_date == "20240510":
    standing_wave = pipeline.medfilter(sp_smoothed.flux, 0.5, 1100/1048576)
    sw = Spectrum1D(standing_wave * u.K,
                    spectral_axis=new_spec_array)

    sp_smoothed = sp_smoothed - sw

# Fitting
init_condition = models.Gaussian1D(-0.2, OH_1665, 0.03*u.MHz)
fit_func = fit_lines(sp_smoothed, init_condition)
sp_fitted = fit_func(sp.frequency)

fit_sp = Spectrum1D(flux=sp_fitted * u.K,
                    spectral_axis=new_spec_coord)
velo_sp = Spectrum1D(flux=sp_fitted * u.K,
                     spectral_axis=new_spec_coord.to(u.km/u.s))


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
# plt.plot(sp.frequency, baseline, label="baseline")
# plt.plot(sp.frequency, sp_baseline_removed.flux, label="baseline removed")
plt.plot(sp.frequency, sp_smoothed.flux, label="smoothed")
# plt.plot(sp.frequency, standing_wave, label="standing wave")
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
velo_array = fit_sp.velocity

plt.figure(figsize=(10,5))
plt.plot(velo_array, sp_smoothed.flux, label="smoothed")
plt.plot(velo_array, fit_sp.flux, label="fitted")


v_min = velo_array[-1].value
v_max = velo_array[0].value
plt.xlim([v_min, v_max])
plt.xlabel("Velocity (km/s)", fontsize=16)
plt.ylabel("Ta (K)", fontsize=16)
plt.grid()
plt.title(select_date, fontsize=20)
plt.legend()
plt.show()
