import copy
import glob
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SpectralCoord
from specutils import Spectrum1D
from specutils.manipulation import FluxConservingResampler
from uwb_tool import pipeline
from uwb_tool import functions


# based_data_path = "./"
based_data_path = "F:/comet_2024/"
comet_name = "12P"
# receiver = "UWB2"
# fixed_freq_limit = [1550, 1750]
receiver = "UWB4"
fixed_freq_limit = [3300, 3450]

obs_date_list = ["20240417", "20240424",
                 "20240510", "20240511", "20240513"]


# prepare for the frequency axis
_, freq_array = functions.get_freq_mask_range(receiver, fixed_freq_limit[0], fixed_freq_limit[1])

c = 3e5 # light velocity
source_on_time = 5 # minutes
obs_freq_array = freq_array
new_spec_array = copy.deepcopy(freq_array) * u.MHz


# read the ta calibrated ta data files
ta_on_files = f"{based_data_path}" + "comet_outputs/{}/{}/{}/tacal/sourceON/Ta_{}_{}_*.npy"
ta_off_files = f"{based_data_path}" + "comet_outputs/{}/{}/{}/tacal/sourceOFF/Ta_{}_{}_*.npy"

ta_on, ta_off = {}, {}
for obs_date in obs_date_list:
    ta_on[obs_date] = None
    ta_off[obs_date] = None

ta_onoff = {}
for obs_date in obs_date_list:
    ta_onoff[obs_date] = None


for obs_date in obs_date_list:
    for ta_on_file in sorted(glob.glob(ta_on_files.format(comet_name, obs_date, receiver, fixed_freq_limit[0], fixed_freq_limit[1]))):
        temp_ta = np.load(ta_on_file)
        temp_ta = temp_ta
        if ta_on[obs_date] is None:
            ta_on[obs_date] = temp_ta
        else:
            ta_on[obs_date] = np.vstack((ta_on[obs_date], temp_ta))
    
    for ta_off_file in sorted(glob.glob(ta_off_files.format(comet_name, obs_date, receiver, fixed_freq_limit[0], fixed_freq_limit[1]))):
        temp_ta = np.load(ta_off_file)
        temp_ta = temp_ta
        if ta_off[obs_date] is None:
            ta_off[obs_date] = temp_ta
        else:
            ta_off[obs_date] = np.vstack((ta_off[obs_date], temp_ta))

# remove the last un-paired on-source observations
for obs_date in obs_date_list:
    ta_on[obs_date] = ta_on[obs_date][0:-5,:]

# remove the glitch on 20240510 of the last 5 minutes in off-source
ta_on["20240510"] = ta_on["20240510"][0:-5,:]
ta_off["20240510"] = ta_off["20240510"][0:-5,:]

# doppler correction for on-source data
for obs_date in obs_date_list:
    ta_onoff[obs_date] = np.zeros_like(ta_on[obs_date])

    # read the ephem file
    ephem_file = f"./comet_ephem/horizons_results_{comet_name}_{obs_date[-4:]}_more.txt"
    comet_velo_interp = pipeline.comet_ephem(ephem_file, tstr2mjd=False)

    onoff_cycle = -1
    for idx in range(ta_on[obs_date].shape[0]):
        functions.prt_info("Doppler correction for sourceON on %s of %d", obs_date, idx)
        if idx % source_on_time == 0:
            onoff_cycle += 1
        relative_time = idx + onoff_cycle*source_on_time
        comet_velo = comet_velo_interp(relative_time)

        rest_freq_array = obs_freq_array / (1.0 - comet_velo/c)
        
        temp_ta = ta_on[obs_date][idx, :] - ta_off[obs_date][idx, :]

        sp = Spectrum1D(flux=temp_ta * u.K,
                        spectral_axis=rest_freq_array * u.MHz,
                        velocity_convention="radio")

        # resampling
        flux_conservation_resample = FluxConservingResampler()
        sp_resampled = flux_conservation_resample(sp, new_spec_array)
        
        ta_onoff[obs_date][idx] = sp_resampled.flux


# doppler correction for off-source data
# for obs_date in obs_date_list:
#     # read the ephem file
#     ephem_file = f"./comet_ephem/horizons_results_{comet_name}_{obs_date[-4:]}_more.txt"
#     comet_velo_interp = pipeline.comet_ephem(ephem_file, tstr2mjd=False)
#     onoff_cycle = 0
#     for idx, temp_ta in enumerate(ta_off[obs_date]):
#         functions.prt_info("Doppler correction for sourceOFF on %s of %d", obs_date, idx)
#         if idx % source_on_time == 0:
#             onoff_cycle += 1
#         comet_velo = comet_velo_interp(idx + onoff_cycle*source_on_time)
#         rest_freq_array = obs_freq_array / (1.0 - comet_velo/c)
        
#         sp = Spectrum1D(flux=temp_ta * u.K,
#                         spectral_axis=rest_freq_array * u.MHz,
#                         velocity_convention="radio")

#         # resampling
#         flux_conservation_resample = FluxConservingResampler()
#         sp_resampled = flux_conservation_resample(sp, new_spec_array)
        
#         ta_off[obs_date][idx] = sp_resampled.flux

# on-source - off-source
# ta_onoff = {}
# for obs_date in obs_date_list:
#     ta_onoff[obs_date] = ta_on[obs_date] - ta_off[obs_date]

# average ta data at different time
for obs_date in obs_date_list:
    ta_onoff[obs_date] = np.average(ta_onoff[obs_date], axis=0)

# save the data
for obs_date in obs_date_list:
    functions.prt_info("Writing doppler corrected data on %s", obs_date)
    
    out_path = f"{based_data_path}/comet_outputs/{comet_name}/{obs_date}/{receiver}/product"
    ta_doppler_file = f"Ta_{fixed_freq_limit[0]}_{fixed_freq_limit[1]}_doppler.npy"
    
    np.save(f"{out_path}/sourceON/{ta_doppler_file}", ta_on[obs_date])
    np.save(f"{out_path}/sourceOFF/{ta_doppler_file}", ta_off[obs_date])
    np.save(f"{out_path}/{ta_doppler_file}", ta_onoff[obs_date])
