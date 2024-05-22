import glob
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.modeling import models
from astropy.coordinates import SpectralCoord
from specutils import Spectrum1D
from specutils.manipulation import box_smooth, FluxConservingResampler
from specutils.fitting.continuum import fit_continuum
from specutils.fitting import fit_lines
# from uwb_tool import pipeline
from uwb_tool import functions


comet_name = "12P"
obs_date_list = ["20240417", "20240424",
                 "20240429", "20240503",
                 "20240510", "20240511", "20240513"]
receiver = "UWB2"
fixed_freq_limit = [1655, 1670]
freq_limit = [1664, 1668]
OH_1667 = 1667.3590 * u.MHz
select_date = obs_date_list[0]


tacal_on_files = "./comet_outputs/{}/{}/tacal/sourceON/Ta_{}_{}_*.npy"
tacal_off_files = "./comet_outputs/{}/{}/tacal/sourceOFF/Ta_{}_{}_*.npy"

ta_on, ta_off = {}, {}
for obs_date in obs_date_list:
    ta_on[obs_date] = None
    ta_off[obs_date] = None

for obs_date in obs_date_list:
    for tacal_on_file in sorted(glob.glob(tacal_on_files.format(comet_name, obs_date, fixed_freq_limit[0], fixed_freq_limit[1]))):
        temp_tacal = np.load(tacal_on_file)
        # print(tacal_on_file)
        # print(temp_tacal.shape)
        if ta_on[obs_date] is None:
            ta_on[obs_date] = temp_tacal
        else:
            ta_on[obs_date] = np.vstack((ta_on[obs_date], temp_tacal))
    
    for tacal_off_file in sorted(glob.glob(tacal_off_files.format(comet_name, obs_date, fixed_freq_limit[0], fixed_freq_limit[1]))):
        temp_tacal = np.load(tacal_off_file)
        if ta_off[obs_date] is None:
            ta_off[obs_date] = temp_tacal
        else:
            ta_off[obs_date] = np.vstack((ta_off[obs_date], temp_tacal))

print(ta_on[select_date].shape,ta_off[select_date].shape)
exit(0)



_, freq_array = functions.get_freq_mask_range(receiver, fixed_freq_limit[0], fixed_freq_limit[1])
freq_mask = (freq_array >= freq_limit[0]) & (freq_array <= freq_limit[1])

ta = ta_onoff[select_date][freq_mask]

freq_axis = SpectralCoord(freq_array[freq_mask] * u.MHz,
                          doppler_convention="radio",
                          doppler_rest=OH_1667)

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
v_lsr = freq_axis.to(u.km / u.s, doppler_convention='radio', doppler_rest=OH_1667)

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