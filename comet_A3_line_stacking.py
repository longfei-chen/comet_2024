import sys
import copy
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.modeling import models
from astropy.coordinates import SpectralCoord
# from astropy.nddata import StdDevUncertainty
from specutils import Spectrum1D
# from specutils import SpectralRegion
from specutils.analysis import gaussian_fwhm
# from specutils.analysis import centroid, snr_derived, line_flux
from specutils.manipulation import box_smooth
from specutils.fitting.continuum import fit_continuum
from specutils.fitting import fit_lines

from uwb_tool import functions
from uwb_tool import linedb
from uwb_tool import molspec
import utils


selected_freq_list = [
[1300,1350], [1350,1450], [1630,1664], [1664,1670], [1705,1750],
[1750,1800], [1975,2060],[2060,2100], [2100,2150], [2235,2325],
[2600,2660], [2700,2750], [2750,2845], [2870,2940], [2960,3000], 
[3050,3120], [3120,3150], [3150,3300], [3320,3450]
]

comet = "A3"
molecule = sys.argv[1]
obsdates = list(utils.comet_property[comet].keys())
data_path = "/media/longfei/c57a45a4-0626-4f7a-bdbb-f0d65a153c9d/N2024_8_selected/"


line_list = linedb.from_splatalogue("./line_db/splatalogue_1.2-3.5GHz.csv")
line_list.update(utils.line_db)
# print(line_list.keys())
# exit(0)

rest_frequency_list = line_list[molecule]

rest_frequency_list = functions.clean_restfreqlist(rest_frequency_list, selected_freq_list, nearest=1.0)
# print(rest_frequency_list)
# exit(0)

comet_data = molspec.molspec(comet=comet, mol=molecule, 
                             restfreqs=rest_frequency_list, obsdates=obsdates,
                            loadpath=data_path)

mole_lines = comet_data.lines

# mole_lines.show_spec()
# mole_lines.show_ss()

for rest_freq, specs in zip(mole_lines.restfreqs , mole_lines.speclist):
    # functions.prt_info("Plotting spectrum for transition %f MHz", rest_freq)
    # spec.show_spec()

    averaged_ta = None
    
    for obs_date, spec_data in zip(specs.obsdates, specs.specs):
        # if obs_date != "20241013" or int(rest_freq) != 1667:
        #     continue

        freq_arr = spec_data[:, 0]
        ta = spec_data[:, 2]

        new_spec_array = freq_arr * u.GHz
        new_spec_coord = SpectralCoord(new_spec_array,
                                    doppler_convention="radio",
                                    doppler_rest=rest_freq*u.MHz)
        
        sp = Spectrum1D(flux=ta * u.K,
                        spectral_axis=new_spec_coord,
                        velocity_convention="radio")

        # Baseline removing
        fitted_continuum = fit_continuum(sp)
        baseline = fitted_continuum(sp.frequency)
        sp_baseline_removed = sp - baseline

        # Smoothing
        sp_smoothed = box_smooth(sp_baseline_removed, width=5)

        if averaged_ta is None:
            averaged_ta = sp_smoothed.flux
        else:
            averaged_ta += sp_smoothed.flux
        # continue

        # Fitting
        init_condition = models.Gaussian1D(-0.2, rest_freq*u.MHz, 0.03*u.MHz)
        fit_func = fit_lines(sp_smoothed, init_condition)
        sp_fitted = fit_func(sp.frequency)
        

        velo_sp = Spectrum1D(flux=sp_fitted * u.K, 
                            spectral_axis=new_spec_coord.to(u.km/u.s))


        rms = functions.get_rms(sp_smoothed.flux, sigma=3)
        peak_ta_idx = np.argmax(np.abs(sp_fitted))
        peak_ta = sp_fitted[peak_ta_idx]


        functions.prt_info("Showing spectrum for %s at transition %f MHz", obs_date, rest_freq)
        # print(f"SNR: {snr_derived(sp_smoothed, SpectralRegion(1.664*u.GHz, 1.668*u.GHz))}")
        print(f"{'rms (mK)':16s}{'peak Ta (mK)':16s}{'SNR':16s}{'FWHM (km/s)':16s}{'center velocity (km/s)'}")
        print(f"{rms.to(u.mK).value:<16.2f}"
              f"{peak_ta.to(u.mK).value:<16.2f}"
              f"{np.abs(peak_ta/rms):<16.2f}"
              f"{gaussian_fwhm(velo_sp).value:<16.2f}"
              f"{velo_sp.velocity[peak_ta_idx].value:<16.2f}")


        integ_flux = integrate.simpson(y=velo_sp.flux, x=velo_sp.velocity)
        gain = utils.get_gain(rest_freq)
        integ_flux_Jy = integ_flux / gain

        # save_data = np.array([velo_sp.velocity.value, sp.flux.value/gain]).T
        # np.save("comet_A3_20241013.npy", save_data)


        velo_reso = velo_sp.velocity[0] - velo_sp.velocity[1]
        integ_flux_std_dev = rms.value * np.sqrt(integ_flux_Jy * velo_reso.value)

        earth_comet_dist = utils.comet_property[comet][obs_date][0] * u.au
        sun_comet_dist = utils.comet_property[comet][obs_date][1] * u.au
        inversion_factor = abs(utils.comet_property[comet][obs_date][-1])

        Q = utils.Q_Bockellee1990(earth_comet_dist, sun_comet_dist, integ_flux_Jy*u.Jy*u.km/u.s, inversion=inversion_factor)
        Q_err = utils.Q_Bockellee1990(earth_comet_dist, sun_comet_dist, integ_flux_std_dev*u.Jy*u.km/u.s, inversion=inversion_factor)

        print(f"Integ_flux (Jy*km/s): {integ_flux_Jy:<10.3e} +- {integ_flux_std_dev:<.3e}")
        print(f"Q_{molecule} (s^-1): {Q:<10.3e} +- {Q_err:.3e}")



        # plot with respect to frequency
        plt.figure(figsize=(10,5))
        plt.plot(sp.frequency, sp_smoothed.flux, label="smoothed")
        plt.plot(sp.frequency, sp_fitted, label="fitted")


        freq_min = sp.frequency[0].value
        freq_max = sp.frequency[-1].value
        plt.xlim([freq_min, freq_max])
        plt.xlabel("Frequency (GHz)", fontsize=16)
        plt.ylabel("Ta (K)", fontsize=16)
        plt.grid()
        plt.title(f"{obs_date} @ {molecule} {rest_freq} MHz", fontsize=20)
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
        plt.title(f"{obs_date} @ {molecule} {rest_freq} MHz", fontsize=20)
        plt.legend()
        plt.show()
    
    # show averaged spectrum
    averaged_ta /= len(specs.specs)
    specs.intensity[:,2] = averaged_ta

    specs.intensity = specs.int_cutoff(vmin=-30, vmax=30)

    rms = functions.get_rms(specs.intensity[:,2])

    functions.prt_info("-"*40)
    functions.prt_info("Showing averaged spectrum at transition %f MHz for all observed dates", rest_freq)
    print(f"averaged spectrum rms: {rms*1e3:.2f} mK")

    specs.show_int(xunit="velo")


mole_lines.stacked_spec = mole_lines.stacking()

rms = functions.get_rms(mole_lines.stacked_spec[:, 1])

functions.prt_info("-"*40)
functions.prt_info("Showing stacked spectrum for all transitions")
print(f"stacked spectrum rms: {rms*1e3:.2f} mK")

mole_lines.show_ss()
