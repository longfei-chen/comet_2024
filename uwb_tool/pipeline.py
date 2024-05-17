#!/usr/bin/python3
#-*- encoding: utf-8 -*-

'''
pipeline for spectral line process for FAST.
'''
import os, gc, glob
import numpy as np
from astropy.time import Time
from astropy.modeling import models, fitting

from scipy.io.idl import readsav
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import scipy.signal as signal

# from specline import constants, molspec
from uwb_tool import functions, plotting
from uwb_tool.ugdopplerfast import ugdopplerfast




def polave_fits(obssetup, src_path, out_path, fig_path=None):
    '''
    average the XX, YY polarization data.

    input: obssetup, src_path, out_path, fig_path=None
        obssetup: configure dict
        src_path: string. raw FAST fits file path.
        out_path: string. polarization averaged file.
    
    return: none
        produce the polarazed power data file.
        output file example: power_0000.npy
    '''
    #check the exists of the source file
    if not functions.is_dir_exists(src_path):
        functions.prt_info("path %s is not exist.", src_path)
        exit(0)
    
    #check and make the output file path
    functions.is_dir_exists(out_path, mkdir=True)
    if fig_path != None:
        functions.is_dir_exists(fig_path, mkdir=True)

    out_file_fmt = "_".join(["power", "{}.npy"])
    
    fits_files = sorted(glob.glob("/".join([src_path, obssetup["fits_files"]])))
    for fi in fits_files:
        functions.prt_info("processing file %s...", fi)
        
        file_idx = int(fi[-9:-5])

        hdu_data = functions.get_fits_data(fi)

        if fig_path != None:
            xx = np.average(hdu_data[:, :, 0], axis=0)
            yy = np.average(hdu_data[:, :, 1], axis=0)

            functions.prt_info("plot polarization...")
            plotting.plot_pol(fig_path, file_idx, xx, yy)
            del xx, yy

        pol_averaged = np.average(hdu_data[:, :, 0:2], axis=2)

        out_file = out_file_fmt.format("%0.4d"%file_idx)
        out_file = out_path + "/" + out_file
        np.save(out_file, pol_averaged)
        
        del pol_averaged, hdu_data
        gc.collect()

    functions.prt_info("Combined polarization finished for %s...", src_path)


    
def power_ave(load_path, prod_path):
    '''
    average the power data of the session.

    input: load_path, prod_path
        load_path: string. polarized power data path 
        prod_path: string. averaged polarized power data path

    return: none
        produce the averaged power data file of the session.
        example: power_ave.npy
    '''
    #check the exists of the data path
    if not functions.is_dir_exists(load_path):
        functions.prt_info("path %s is not exist.", load_path)
        exit(0)

    #check and make the out path dir
    functions.is_dir_exists(prod_path, mkdir=True)

    sum_power = 0
    tot_timebin = 0
    for fi in sorted(glob.glob("/".join([load_path, "power_*.npy"]))):
        functions.prt_info("average power data for file %s...", fi)

        power = np.load(fi)
        tot_timebin += len(power)

        sum_power += np.sum(power, axis=0)

    power_ave = sum_power / tot_timebin

    out_file = "_".join(["power", "ave.npy"])

    functions.prt_info("saving file %s to %s...", out_file, prod_path)
    np.save(prod_path + "/" + out_file, power_ave)
    
    functions.prt_info("Finished the averaged power data for %s...", load_path)


def ta_ave(obssetup, data_path, prod_path, freqlimit=None):
    '''
    average the Ta data of the session.

    input: obssetup, data_path, prod_path, freqlimit
        obssetup: configure dict
        data_path: string. calibrated Ta data path
        prod_path: string. averaged Ta data path
        freqlimit: list. frequency range, example: [1419, 1421]

    return: none
        produce the averaged Ta data file of the session.
        output file example: Ta_1419_1421_ave.npy
    '''
    #check the exists of the data path
    if not functions.is_dir_exists(data_path):
        functions.prt_info("path %s is not exist.", data_path)
        exit(0)

    #check and make the out path dir
    functions.is_dir_exists(prod_path, mkdir=True)

    if freqlimit is None:
        freqlimit = functions.UWB_receiver[obssetup["receiver"]]["recommended_range"]

    if (freqlimit[0] < functions.UWB_receiver[obssetup["receiver"]]["total_range"][0]) \
        or (freqlimit[1] > functions.UWB_receiver[obssetup["receiver"]]["total_range"][1]):
        functions.prt_info("freqlimit is out of the receiver frequency range.")
        exit(0)

    sum_ta = 0
    tot_record = 0
    for fi in sorted(glob.glob("/".join([data_path, f"Ta_{freqlimit[0]}_{freqlimit[1]}_*.npy"]))):
        functions.prt_info("average Ta data for file %s...", fi)

        ta = np.load(fi)
        tot_record += 1
        
        sum_ta += ta

    ta_ave = sum_ta / tot_record

    out_file = "_".join(["Ta", str(freqlimit[0]), str(freqlimit[1]), "ave.npy"])

    functions.prt_info("saving file %s to %s...", out_file, prod_path)
    np.save(prod_path + "/" + out_file, ta_ave)
    
    functions.prt_info("Finished the averaged Ta data for %s...", data_path)
    
    
def save_molspec(molspec, filename):
    np.save(molspec, filename)

def load_molspec(filename):
    return np.load(filename)
    
    
def medfilter(data_arr, window, resolution=500/1048576):
    '''
    median filter for the spectral.

    input: data_arr, window, resolution=500/1048576
        data_arr: list. Ta
        window: filter window in MHz
        resolution: float in MHz / per channel

    return: fitting array.
    '''
    box_width = int(window / resolution)
    if box_width % 2 == 0:
        box_width += 1

    baseline = signal.medfilt(data_arr, kernel_size=box_width)

    return baseline

#def polyfit(x_arr, y_arr, order=3):
#    coeff = np.polyfit(x_arr, y_arr, order)
#    fitfunc = np.poly1d(coeff)
#    
#    return fitfunc(x_arr)

def polyfit(x_arr, y_arr, order=3):
    poly_model = models.Polynomial1D(degree=order)
    fitter = fitting.LevMarLSQFitter()
    poly_fit = fitter(poly_model, x_arr, y_arr)
    
    return poly_fit(x_arr)

def vlsr_corr(target_ra, target_dec, mjd):
    '''
    velocity correction with respect to the target.

    input: target_ra, target_dec, mjd
    target_ra: float in degree.
    target_dec: float in degree.
    mjd: float.

    return: float in km/s
        corrected velocity with respect to the target.
    '''
    jd = mjd + 2400000.5

    vlsr = ugdopplerfast([target_ra], [target_dec], [jd])

    return vlsr[0]

def rest_corr(vlsr, restfreq):
    '''
    doppler correction of the rest frequency.

    input: vlsr, restfreq
    vlsr: float in km/s. corrected velocity of the observer.
    restfreq: float in MHz. rest frequency of the molecular line.

    return: float in MHz
        corrected doppler of the molecular line.
    '''
    doppler = vlsr * restfreq / 2.99792458e5

    return doppler

def freq_shift(freq_arr, doppler, resolution=500/1048576):
    '''
    correction of the original frequency sequency.

    input: freq_arr, doppler, resolution=500/1048576
    freq_arr: original frequency sequency.
    doppler: corrected doppler of the molecular line.

    return: array
    corrected original frequency sequency.
    '''
    freq_arr += doppler

    channel_offset = int(doppler / resolution)

    return functions.shift_arr(freq_arr, channel_offset)

def ta_shift(ta_arr, doppler, resolution=500/1048576):
    '''
    correction of the original Ta sequency.

    input: ta_arr, doppler, resolution=500/1048576
    ta_arr: original Ta sequency.
    doppler: corrected doppler of the molecular line.

    return: array
    corrected original Ta sequency.
    '''
    channel_offset = int(doppler / resolution)

    return functions.shift_arr(ta_arr, channel_offset)


def freq_to_velo(freq_arr, restfreq):
    '''
    convertion from frequency domain to velocity domain.

    input: freq_arr, restfreq
    freq_arr: corrected frequency.
    restfreq: rest molecular frequency.

    return: array.
    frequency array to velocity array.
    '''
    velo = 2.99792458e5 * (restfreq - freq_arr) / restfreq

    return velo

def ta2flux(ta_arr, eta=0.53):
    '''
    convertion from Ta to flux.

    input: ta_arr, eta=0.53
    ta_arr: antenna temperature.

    return: array in Jy.
    '''
    functions.prt_info("flux calibration...")

    gain = eta * 25.2

    return ta_arr / gain


def merge_figs(obssetup, fig_path, prod_path, freqlimit=None, fig_pattern_list=["*_dyn.pdf", "*_Ta.pdf", "*_timebin.pdf"]):
    if not functions.is_dir_exists(fig_path):
        functions.prt_info("path %s is not exist.", fig_path)
        exit(0)
    
    functions.is_dir_exists(prod_path, mkdir=True)
    
    if freqlimit is None:
        freqlimit = functions.UWB_receiver[obssetup["receiver"]]["recommended_range"]

    if (freqlimit[0] < functions.UWB_receiver[obssetup["receiver"]]["total_range"][0]) \
        or (freqlimit[1] > functions.UWB_receiver[obssetup["receiver"]]["total_range"][1]):
        functions.prt_info("freqlimit is out of the receiver frequency range.")
        exit(0)

    for figs in fig_pattern_list:
        pdf_figs = fig_path + "/" + "_".join([str(freqlimit[0]), str(freqlimit[1]), figs])
        out_pdf = prod_path + "/" + "_".join([str(freqlimit[0]), str(freqlimit[1]), figs[2:]])
        
        merge_command = "pdfunite {} {}".format(pdf_figs, out_pdf)
        
        functions.prt_info("merging figures for %s...", pdf_figs)
        try:
            os.system(merge_command)
        except:
            functions.prt_info("The shell commend pdfunite is not installed in this system.")
            functions.prt_info("pdf files are not merged.")
            exit(0)


def comet_split_source_onoff(obssetup, load_path, out_path):
    """
    output file name example: power_source_on_0001.npy
    """
    #check the exists of the load file
    if not functions.is_dir_exists(load_path):
        functions.prt_info("path %s is not exist.", load_path)
        exit(0)
    
    #check and make the output data path
    source_on_path = out_path + "/" + "sourceON"
    source_off_path = out_path + "/" + "sourceOFF"

    functions.is_dir_exists(source_on_path, mkdir=True)
    functions.is_dir_exists(source_off_path, mkdir=True)
    
    source_on_file_fmt = "_".join(["power", "source_on", "{}.npy"])
    source_off_file_fmt = "_".join(["power", "source_off", "{}.npy"])

    sampling_time = obssetup["sampling_time"]
    source_on_time = obssetup["source_on"]
    source_switch_time = obssetup["source_switch"]

    split_onoff_idx = int(source_on_time/sampling_time)
    split_switch_idx = int((source_on_time + source_switch_time)/sampling_time)
    
    source_onoff_switch = -1

    idx = 0
    power = None
    power_files = sorted(glob.glob("/".join([load_path, "power_*.npy"])))
    the_last_file = False
    for fi_idx,fi in enumerate(power_files):
        if fi_idx+1 == len(power_files):
            the_last_file = True

        functions.prt_info("reading power file %s...", fi)

        if power is None:
            power = np.load(fi)
        else:
            foo = np.load(fi)

            power = np.concatenate((power, foo), axis=0)
            del foo

        #read the data stream until some cycle of source onoff is done
        if (not the_last_file) and (len(power)*sampling_time < source_on_time + source_switch_time):
            continue
        
        while the_last_file or (len(power)*sampling_time > source_on_time + source_switch_time):
            if the_last_file:
                the_last_file = False

            functions.prt_info("splitting source ON and OFF power data...")
            functions.prt_info("writing file to path %s/{sourceON,sourceOFF}...", out_path)

            source_data = power[0:split_onoff_idx, :]
            power = power[split_switch_idx:, :]

            idx += 1
            source_on_file = source_on_file_fmt.format("%0.4d"%idx)
            source_on_file = source_on_path + "/" + source_on_file

            source_off_file = source_off_file_fmt.format("%0.4d"%idx)
            source_off_file = source_off_path + "/" + source_off_file
            
            source_onoff_switch += 1
            if source_onoff_switch % 2 == 0:
                np.save(source_on_file, source_data)
            else:
                np.save(source_off_file, source_data)

            del source_data
            gc.collect()
        # end while

def comet_cal_power_data(obssetup, load_path, data_path, freqlimit=None, cycle=3, fig_path=None):
    """
    output file name example: Ta_1665_1667_0000.npy
    """
    #check the exists of the load file
    if not functions.is_dir_exists(load_path):
        functions.prt_info("path %s is not exist.", load_path)
        exit(0)
    
    #check and make the output file path
    functions.is_dir_exists(data_path, mkdir=True)
    if fig_path != None:
        functions.is_dir_exists(fig_path, mkdir=True)
    
    if freqlimit is None:
        freqlimit = functions.UWB_receiver[obssetup["receiver"]]["recommended_range"]

    if (freqlimit[0] < functions.UWB_receiver[obssetup["receiver"]]["total_range"][0]) \
        or (freqlimit[1] > functions.UWB_receiver[obssetup["receiver"]]["total_range"][1]):
        functions.prt_info("freqlimit is out of the receiver frequency range.")
        exit(0)

    out_file_fmt = "_".join(["Ta", str(freqlimit[0]), str(freqlimit[1]), "{}.npy"])

    sampling_time = obssetup["sampling_time"]
    period_on = obssetup["period_on"]
    period_off = obssetup["period_off"]
    period = period_on + period_off

    onoff_steps = int(period/sampling_time)
    sequence = np.arange(onoff_steps*cycle)

    cal_on_mask = sequence % onoff_steps < int(period_on/sampling_time)
    cal_off_mask = ~cal_on_mask

    freq_mask, freq_range = functions.get_freq_mask_range(obssetup["receiver"], xmin=freqlimit[0], xmax=freqlimit[1])
    
    tcal_interp1d = functions.get_tcal(obssetup["tcal_file"], pol_ave=True)
    Tcal = tcal_interp1d(freq_range)

    idx = 0
    power = None
    power_files = sorted(glob.glob("/".join([load_path, "power_*.npy"])))
    for fi in power_files:
        functions.prt_info("reading power file %s...", fi)

        power = np.load(fi)
        
        power = power[:, freq_mask]

        #read the source ON or OFF power data to calibrate
        while power.shape[0] != 0:
            idx += 1
            functions.prt_info("calibration for output file %s...", out_file_fmt.format("%0.4d"%idx))

            #split cal on and cal off
            image = power[0:onoff_steps*cycle, :]
            power = power[onoff_steps*cycle:, :]

            if power.shape[0] == 0:
                power_calon = image[cal_on_mask[0:image.shape[0]], :]
                power_caloff = image[cal_off_mask[0:image.shape[0]], :]
            else:
                power_calon = image[cal_on_mask, :]
                power_caloff = image[cal_off_mask, :]

            #average cal on and cal off at time domain
            power_calon = np.average(power_calon, axis=0)
            power_caloff = np.average(power_caloff, axis=0)
            power_calres = power_calon - power_caloff
            
            # if np.average(power_calres) < 0:
            #     functions.prt_info("On-OFF less than 0 at %d sec.", (idx-1)*period*cycle)
            #     functions.prt_info("this file is canceled.")
            #     continue

            Ta = Tcal * power_caloff / np.average(power_calres)
            #Ta = Tcal * power_caloff / power_calres

            #save the calibrated Ta file
            out_file = out_file_fmt.format("%0.4d"%idx)
            out_file = data_path + "/" + out_file
            np.save(out_file, Ta)

            #plot some figures
            if fig_path == None:
                continue
            
            #plot dynamic spectrum
            functions.prt_info("plotting dynamic spectrum...")
            xyrange = [freqlimit[0], freqlimit[1]]
            xyrange += [(idx-1)*period*cycle, idx*period*cycle]
            plotting.plot_dyn(fig_path, freqlimit, idx, xyrange, image)

            #plot timebin
            functions.prt_info("plotting timebin...")
            x = sequence*sampling_time + (idx-1)*period*cycle
            y = np.average(image, axis=1)
            x = x[0:image.shape[0]]
            plotting.plot_timebin(fig_path, freqlimit, idx, period*cycle, x, y)

            #plot Ta
            functions.prt_info("plotting Ta...")
            plotting.plot_Ta(fig_path, freqlimit, idx, freq_range, Ta)
            
            #plot cal
            functions.prt_info("plotting cal...")
            plotting.plot_Cal(fig_path, freqlimit, idx, freq_range, power_calon, "calon")
            plotting.plot_Cal(fig_path, freqlimit, idx, freq_range, power_caloff, "caloff")
            plotting.plot_Cal(fig_path, freqlimit, idx, freq_range, power_calres, "calres")
            
            del power_calon, power_caloff, power_calres, image, Ta, xyrange, x, y
        del power
        gc.collect()
        # end while
    functions.prt_info("Calibration finished for %s...", load_path)

def comet_ON_minus_OFF(obssetup, prod_path, freqlimit=None):
    if not functions.is_dir_exists(prod_path):
        functions.prt_info("path %s is not exits.", prod_path)
        exit(0)
    
    if freqlimit is None:
        freqlimit = functions.UWB_receiver[obssetup["receiver"]]["recommended_range"]

    if (freqlimit[0] < functions.UWB_receiver[obssetup["receiver"]]["total_range"][0]) \
        or (freqlimit[1] > functions.UWB_receiver[obssetup["receiver"]]["total_range"][1]):
        functions.prt_info("freqlimit is out of the receiver frequency range.")
        exit(0)

    sourceON_files = sorted(glob.glob("/".join([prod_path+"/sourceON", f"Ta_{freqlimit[0]}_{freqlimit[1]}_ave.npy"])))
    sourceOFF_files = sorted(glob.glob("/".join([prod_path+"/sourceOFF", f"Ta_{freqlimit[0]}_{freqlimit[1]}_ave.npy"])))

    file_num = min(len(sourceON_files), len(sourceOFF_files))
    for on_file,off_file in zip(sourceON_files[0:file_num], sourceOFF_files[0:file_num]):
        on_file_name = os.path.basename(on_file)
        off_file_name = os.path.basename(off_file)
        if on_file_name != off_file_name:
            functions.prt_info("source ON and OFF file has not the same filename.")
            functions.prt_info("ON file: %s", on_file_name)
            functions.prt_info("OFF file: %s", off_file_name)
            exit(0)
            
        on_data = np.load(on_file)
        off_data = np.load(off_file)
        onoff = on_data - off_data
        
        np.save(prod_path + "/" + on_file_name, onoff)
        functions.prt_info("Baseline substraction for %s...", on_file_name)


def comet_velo_ephem(ephem_file, tstr2mjd=True):
    mjd_arr = []
    velo_arr = []
    
    data_stream = False
    with open(ephem_file, "r", encoding='utf-8') as fp:
        for line in fp.readlines():
            if "$$SOE" in line:
                data_stream = True
                continue
            if "$$EOE" in line:
                break
                
            if data_stream:
                line = line.strip().split()
                
                if tstr2mjd:
                    tstr = " ".join([line[0], line[1]])
                    t = Time(tstr, format="iso", scale="utc")
                    mjd_arr.append(t.mjd)
                else:
                    mjd_arr.append(float(line[0])-2400000.5)
                
                velo_arr.append(float(line[-1]))
            
    return interp1d(mjd_arr, velo_arr, kind="linear")
    
