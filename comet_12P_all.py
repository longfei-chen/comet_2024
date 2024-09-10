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

mod = "onoff"
mol = "CH3CHOHCH2OH"

# catlog = linedb.get_linedb()
catlog = {
"OD": [1190.7741, 1191.1047],
"18OH": [1584.274, 1637.564, 1639.503, 1692.795],
"17OH": [1322.4597, 1455.7225, 1624.5096, 1626.1564, 1656.5422, 1902.0885, 1940.2756, 2008.3517, 2027.319, 2102.7856, 2117.7941],
"NiO": [1924.1406],
"OH": [1612.2310, 1665.4018, 1667.3590, 1720.5300], 
"HI": [1420.4058],
"CH": [3263.794, 3335.481, 3349.193],
"13CH3OH": [794.7061, 2384.0513],
"CH3OCHO": [1610.2493, 1610.9063], 
"HC5N": [2662.6641, 2662.8795],
"HC9N": [2905.1827],
"HC11N": [1690.6293, 2028.7551, 2705.0066, 2705.0066, 3043.1323, 3043.1323],
"CH3CHOHCH2OH": [3349.7184], 
"c-C3H": [3447.7142, 3447.8425, 3447.5665, 3447.6246], 
"H2SO4": [3350.2291],
}
rest_freq_list = catlog[mol]
rest_freq_list = functions.clean_restfreqlist(rest_freq_list, [[500, 3450]])

if len(rest_freq_list) == 0:
    functions.prt_info("No transition lines for %s in the selected frequency range.", mol)
    exit()

load_path = "F:/comet_2024/"
comet_name = "12P"

obs_date_list = ["20240417", "20240424",
                 "20240510", "20240511", "20240513"]

mymolspec = molspec.molspec(obssetup={}, 
                            name=mol, 
                            beam="M01", 
                            restfreqs=rest_freq_list, 
                            mjds=[None]*len(obs_date_list), 
                            obsdates=obs_date_list, 
                            lines=None, 
                            loadpath=None)

mymolspec.lines = molspec.lines()
mymolspec.lines.restfreqs = rest_freq_list
mymolspec.lines.speclist = []

nlines = len(rest_freq_list)
vlsr = [0]*nlines
for iline, eachline in enumerate(mymolspec.restfreqs):
    functions.prt_info("Processing for transition %f MHz for %s (%d/%d)", eachline, mymolspec.name, iline+1, nlines)

    spec = molspec.obsspec()
    spec.obsdates = mymolspec.obsdates
    spec.mjds = mymolspec.mjds
    spec.specs = []

    receiver, fixed_freq_limit = utils.get_band_range(eachline)
    width = 1 # MHz
    freq_limit = [eachline - width, eachline + width]
    _, freq_array = functions.get_freq_mask_range(receiver, fixed_freq_limit[0], fixed_freq_limit[1])
    freq_mask = (freq_array >= freq_limit[0]) & (freq_array <= freq_limit[1])

    new_spec_array = copy.deepcopy(freq_array[freq_mask]) * u.MHz
    new_spec_coord = SpectralCoord(new_spec_array,
                                doppler_convention="radio",
                                doppler_rest=eachline*u.MHz)

    ta_on_file = f"{load_path}" + "comet_outputs/{}/{}/{}/product/sourceON/Ta_{}_{}_doppler.npy"
    ta_off_file = f"{load_path}" + "comet_outputs/{}/{}/{}/product/sourceOFF/Ta_{}_{}_doppler.npy"
    ta_onoff_file = f"{load_path}" + "comet_outputs/{}/{}/{}/product/Ta_{}_{}_doppler.npy"

    for eachobs in spec.obsdates:
        if mod == "on":
            ta_file = ta_on_file.format(comet_name, eachobs, receiver, fixed_freq_limit[0], fixed_freq_limit[1])
        if mod == "off":
            ta_file = ta_off_file.format(comet_name, eachobs, receiver, fixed_freq_limit[0], fixed_freq_limit[1])
        if mod == "onoff":
            ta_file = ta_onoff_file.format(comet_name, eachobs, receiver, fixed_freq_limit[0], fixed_freq_limit[1])
        ta = np.load(ta_file)
        if mod in ["on", "off"]:
            ta = ta.mean(axis=0)
        ta = ta[freq_mask]

        functions.prt_info("Initialization data from file %s...", ta_file)

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
        if eachobs == "20240510":
            standing_wave = pipeline.medfilter(sp_smoothed.flux, 0.5, 1100/1048576)
            sw = Spectrum1D(standing_wave * u.K,
                            spectral_axis=new_spec_coord,
                            velocity_convention="radio")

            sp_smoothed = sp_smoothed - sw

        # # Fitting
        # init_condition = models.Gaussian1D(-0.2, eachline, 0.03*u.MHz)
        # fit_func = fit_lines(sp_smoothed, init_condition)
        # sp_fitted = fit_func(sp.frequency)

        # fit_sp = Spectrum1D(flux=sp_fitted * u.K,
        #                     spectral_axis=new_spec_coord)

        xyz = np.column_stack((sp_smoothed.frequency.to(u.MHz).value, new_spec_coord.to(u.km/u.s).value, sp_smoothed.flux.value))
        spec.specs.append(xyz)
    spec.intensity = spec.ave_spec()
    mymolspec.lines.speclist.append(spec)



lines = mymolspec.lines
spec_list = lines.speclist

for i,spec in enumerate(spec_list):
    print(f"\r{i+1}", end="")
    # spec.show_spec(xunit="freq", vpos=lines.restfreqs[i], title=f"{mol} @{lines.restfreqs[i]} MHz")
    spec.show_spec(xunit="velo", vpos=vlsr[i], title=f"{mol} @{lines.restfreqs[i]} MHz")


##for each rest frequency, average their spectrum for different mjds
# for spec in spec_list:
#    spec.specs = spec.rfi_remove(sigma=5)
#    spec.intensity = spec.ave_spec()
#    spec.intensity = spec.int_cutoff(vmin=-40, vmax=40)


#show intensity profile for each line
# lines.show_spec(xunit="freq")

# lines.stacked_spec = lines.stacking()
# lines.show_ss(title="stacked spectrum for every line")
