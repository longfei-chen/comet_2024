import sys
import copy
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.modeling import models
from astropy.nddata import StdDevUncertainty
from astropy.coordinates import SpectralCoord
from specutils import Spectrum1D, SpectralRegion
from specutils.analysis import centroid, gaussian_fwhm, snr_derived, line_flux
from specutils.manipulation import box_smooth
from specutils.fitting.continuum import fit_continuum
from specutils.fitting import fit_lines

from uwb_tool import pipeline
from uwb_tool import functions
import utils

if len(sys.argv) != 3:
    print(f"Usage: python {sys.argv[0]} OH1665 20240417")
    exit()

line_dict = utils.line_dict

obs_date_list = ["20241003", "20241004", "20241013", "20241014", "20241016", "20241017"]


mol_str = sys.argv[1]
if mol_str not in line_dict:
    print(f"Available lines: {line_dict.keys()}")
    exit()

mol_line = line_dict[mol_str]
width = 1 # MHz
freq_limit = [np.floor(mol_line.value - width), np.floor(mol_line.value + width)]

select_date = str(sys.argv[2])
if select_date not in obs_date_list:
    print(f"Available date: {obs_date_list}")
    exit()


receiver, fixed_freq_limit = utils.get_band_range(mol_line.value)


based_data_path = "/media/longfei/c57a45a4-0626-4f7a-bdbb-f0d65a153c9d/N2024_8_outputs/"
comet_name = "A3"

# prepare for the frequency axis
_, freq_array = functions.get_freq_mask_range(receiver, fixed_freq_limit[0], fixed_freq_limit[1])
freq_mask = (freq_array >= freq_limit[0]) & (freq_array <= freq_limit[1])

new_spec_array = copy.deepcopy(freq_array[freq_mask]) * u.MHz
new_spec_coord = SpectralCoord(new_spec_array,
                               doppler_convention="radio",
                               doppler_rest=mol_line)

# read the doppler corrected ta data files
ta_on_files = f"{based_data_path}" + "{}/{}/{}/product/sourceON/Ta_{}_{}_doppler.npy"
ta_off_files = f"{based_data_path}" + "{}/{}/{}/product/sourceOFF/Ta_{}_{}_doppler.npy"
ta_onoff_files = f"{based_data_path}" + "{}/{}/{}/product/Ta_{}_{}_doppler.npy"

ta_on, ta_off, ta_onoff = {}, {}, {}
ta_onoff_ave = None
for obs_date in obs_date_list:
    ta_on[obs_date] = np.load(ta_on_files.format(comet_name, obs_date, receiver, fixed_freq_limit[0], fixed_freq_limit[1]))
    ta_off[obs_date] = np.load(ta_off_files.format(comet_name, obs_date, receiver, fixed_freq_limit[0], fixed_freq_limit[1]))
    ta_onoff[obs_date] = np.load(ta_onoff_files.format(comet_name, obs_date, receiver, fixed_freq_limit[0], fixed_freq_limit[1]))

    if ta_onoff_ave is None:
        ta_onoff_ave = ta_onoff[obs_date]
    else:
        ta_onoff_ave = np.vstack((ta_onoff_ave, ta_onoff[obs_date]))

ta = ta_onoff[select_date][freq_mask]

# ta_onoff_ave = ta_onoff_ave.mean(axis=0)
# ta = ta_onoff_ave[freq_mask]


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
# if select_date == "20240510":
#     standing_wave = pipeline.medfilter(sp_smoothed.flux, 0.5, 1100/1048576)
#     sw = Spectrum1D(standing_wave * u.K,
#                     spectral_axis=new_spec_array)

#     sp_smoothed = sp_smoothed - sw

# Fitting
init_condition = models.Gaussian1D(-0.2, mol_line, 0.03*u.MHz)
fit_func = fit_lines(sp_smoothed, init_condition)
sp_fitted = fit_func(sp.frequency)


velo_sp = Spectrum1D(flux=sp_fitted * u.K, 
                     spectral_axis=new_spec_coord.to(u.km/u.s))


rms = functions.get_rms(sp_smoothed.flux)
peak_ta_idx = np.argmax(np.abs(sp_fitted))
peak_ta = sp_fitted[peak_ta_idx]


# print(f"SNR: {snr_derived(sp_smoothed, SpectralRegion(1.664*u.GHz, 1.668*u.GHz))}")
print(f"{'rms (mK)':16s}{'peak Ta (mK)':16s}{'SNR':16s}{'FWHM (km/s)':16s}{'center velocity (km/s)'}")
print(f"{rms.to(u.mK).value:<16.2f}"
      f"{peak_ta.to(u.mK).value:<16.2f}"
      f"{np.abs(peak_ta/rms):<16.2f}"
      f"{gaussian_fwhm(velo_sp).value:<16.2f}"
      f"{velo_sp.velocity[peak_ta_idx].value:<16.2f}")

# production rate
# integ_flux = line_flux(velo_sp)
# print(integ_flux)

integ_flux = integrate.simpson(velo_sp.flux, velo_sp.velocity)
gain = utils.get_gain(mol_line.value)
integ_flux_Jy = integ_flux / gain
# Q1 = utils.Q_Drahus2010(integ_flux*u.K*u.km/u.s, 50*u.K, 1.0*u.km/u.s, 1.6*u.au, mol_line, utils.line_int[mol_str])
Q2 = utils.Q_Bockelle1990(0.5*u.au, integ_flux_Jy*u.Jy*u.km/u.s, inversion=0.243)
# Q3 = utils.Q_Haser("90000224", "H2O", "^OH$", 50*u.K, 1.0*u.km/u.s, mol_line, 
#                    "2024-04-17 14:53:00.000", integ_flux*u.K*u.km/u.s, 1e28/u.s)
print(f"{'Integ_flux (K*km/s)':20s}Q_{mol_str} (s^-1)")
print(f"{integ_flux:<20.3e}"
    #   f"{Q1:<10.3e}"
      f"{Q2:<10.3e}")


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
plt.title(f"{select_date} @ {mol_str} {mol_line}", fontsize=20)
plt.legend()
plt.show()


# plot with respect to velocity
velo_array = velo_sp.velocity

plt.figure(figsize=(10,5))
plt.plot(velo_array, sp_smoothed.flux, label="smoothed")
plt.plot(velo_array, velo_sp.flux, label="fitted")


v_min = velo_array[-1].value
v_max = velo_array[0].value
plt.xlim([v_min, v_max])
plt.xlabel("Velocity (km/s)", fontsize=16)
plt.ylabel("Ta (K)", fontsize=16)
plt.grid()
plt.title(f"{select_date} @ {mol_str} {mol_line}", fontsize=20)
plt.legend()
plt.show()
