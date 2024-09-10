#!/usr/bin/python3
#-*- encoding: utf-8 -*-

'''
molecule spectral class
'''
import os
import numpy as np
import matplotlib.pyplot as plt

from uwb_tool import pipeline
from uwb_tool import functions
from uwb_tool import plotting


class obsspec():
    '''
    obsspec class to store the observational data for a list of mjds.

    attributes: obsdates, mjds, specs, intensity
    obsdates: list of string. observational date.
        mjds: list of float. observational mjds.
       specs: list of 2d array of the spectral data.
              column one: frequency
              column two: velocity
              column three: Ta
              len(obsdate) == len(mjds) == len(specs)
    intensity: 2d array of the averaged spectral data.
               column one: frequency
               column two: velocity
               column three: Ta
    
    methods:
        ave_spec()
        int_cutoff(vmin=-15, vmax=15, vstep=0.101)
        rfi_remove(sigma=6)
        show_spec(mjd_idx=[1,2,3], xunit="velo | freq")
        show_int(xunit="velo | freq")
    '''
    def __init__(self, obsdates=[], mjds=[], specs=[]):
        self.obsdates = obsdates
        self.mjds = mjds
        self.specs = specs
        
        self.intensity = self.ave_spec()
        
        self.__init_check()
        
    
    def __init_check(self):
        if len(self.obsdates) != len(self.mjds) and len(self.mjds) != len(self.specs):
            functions.prt_info("Length of obsdates, mjds and specs should be equal.")
            return

    def __attr_check(self, primary=False):
        if len(self.specs) == 0:
            functions.prt_info("No data is initialized for obsspec.specs")
            return
        
        for spec in self.specs:
            if type(spec) is not np.ndarray:
                functions.prt_info("Invailed type in the obsspec.speclist, should be <class 'numpy.ndarray'>.")
                return

        #only check the main attribute
        if primary:
            return
            
        if self.intensity is None:
            functions.prt_info("Spectral average is not performed.")
            functions.prt_info("Try: self.ave_spec()")
            return

    def ave_spec(self):
        '''
        average all the mjds spectral.

        input: none

        return: 2d array of the spectral data.
        '''
        self.__attr_check(primary=True)

        if len(self.mjds) == 0:
            return None

        intensity = 0
        N = len(self.mjds)
        
        for spec in self.specs:
            intensity += spec

        intensity /= N

        return intensity
    
    def int_cutoff(self, vmin=-15, vmax=15, vstep=0.101):
        '''
        intensity truncation to a specified velocity range.
        
        input: vmin=-15, vmax=15, vstep=0.101
        
        return: 2d array of the spectral data.
        '''
        self.__attr_check()
        
        velo_arr = self.intensity[:, 1]
        velo_mask = (velo_arr >= vmin) & (velo_arr <= vmax)
        
        profile = self.intensity[velo_mask]
        
        freq_arr = profile[:, 0]
        velo_arr = profile[:, 1]
        ta_arr = profile[:, 2]
        
        velo_interp = np.arange(vmin, vmax, vstep)
        
        freq_interp = np.interp(velo_interp, velo_arr[-1::-1], freq_arr[-1::-1], left=np.nan, right=np.nan)
        ta_interp = np.interp(velo_interp, velo_arr[-1::-1], ta_arr[-1::-1], left=np.nan, right=np.nan)
        
        return np.column_stack((freq_interp, velo_interp, ta_interp))
        
    def rfi_remove(self, sigma=6):
        self.__attr_check()
        
        specs_new = []
        for spec in self.specs:
            rms = functions.get_rms(spec[:,2])
            
            spec[spec[:,2] > sigma*rms,2] = np.nan
            spec[spec[:,2] < -sigma*rms,2] = np.nan
            
            specs_new.append(spec)
        
        return specs_new
        
    
    def show_spec(self, mjd_idx=None, vpos=None, xunit="velo", overlay=False, tight=False, title=None):
        '''
        plot the spectral line(s) for selected mjd(s).

        input: mjd_idx=None, xunit="velo", overlay=False
        mjd_idx: selected mjd to plot.
                 example: mjd_idx = None, 1, or [1,2,3]
        xunit: type of the x-axis of the spectral. velocity or frequency
              example: xunit = "velo" or "freq"
        overlay: boolean
        tight: boolean

        return: none
        '''
        self.__attr_check(primary=True)

        #check the vaildation of the inputs.
        if mjd_idx is None:
            mjd_idx = range(len(self.mjds))
        elif type(mjd_idx) == int:
            mjd_idx = [mjd_idx]
        elif type(mjd_idx) != list or type(mjd_idx) != type(np.array([])):
            functions.prt_info("Invalid type for mjd_idx, must be list or np.array")
            return
        elif len(mjd_idx) == 0:
            mjd_idx = range(len(self.mjds))
        
        for idx in mjd_idx:
            if type(idx) != int:
                functions.prt_info("Invalid type for index, must be int.")
                return
            if idx < 0 or idx >= len(self.mjds):
                functions.prt_info("Input index for mjds out of range. It should be less than %d", len(self.mjds))
                return

        if xunit == "freq":
            sel_col = 0
        elif xunit == "velo":
            sel_col = 1
        else:
            functions.prt_info("Unknown value for xunit: %s", xunit)
            functions.prt_info("Only accept 'freq' or 'velo'.")
            return

        #plot grid figures
        x_arr_list, y_arr_list, label_list = [], [], []

        for idx in mjd_idx:
            spec = self.specs[idx]
            x_arr_list.append(spec[:, sel_col])
            y_arr_list.append(spec[:, 2])
            label_list.append(self.obsdates[idx])

        plotting.plot_spec(x_arr_list, y_arr_list, label_list, vpos=vpos, xunit=xunit, overlay=overlay, tight=tight, title=title)
        

    def show_int(self, xunit="velo", vpos=None, title=None):
        '''
        plot the averaged spectral of all the spectral lines.

        input: xunit="velo"
        xunit: type of the x-axis of the spectral. velocity or frequency
              example: xunit = "velo" or "freq"

        return: none
        '''
        self.__attr_check()

        if xunit == "freq":
            sel_col = 0
        elif xunit == "velo":
            sel_col = 1
        else:
            functions.prt_info("Unknown value for xunit: %s", xunit)
            functions.prt_info("Only accept 'freq' or 'velo'.")
            return

        x_arr = self.intensity[:, sel_col]
        y_arr = self.intensity[:, 2]
        label = "average intensity"
        
        plotting.plot_spec([x_arr], [y_arr], [label], vpos=vpos, xunit=xunit, title=title)


class lines():
    '''
    lines class to store the molecular line transitions.

    attributes: restfreqs, speclist, stacked_spec
    restfreqs: list of rest frequencies.
    speclist: list of obsspec class.
              len(restfreqs) == len(speclist)
    stacked_spec: 2d array of stacked obsspec:
                  column one: velocity
                  column two: intensity of Ta
    mf: 2d array of the matched filter
        column one: velocity
        column two: SNR
    
    methods:
        stacking(by="natural | rms")
        matched_filter(model)
        show_spec()
        show_ss()
        show_mf()
    '''
    def __init__(self, restfreqs=[], speclist=[]):
        self.restfreqs = restfreqs
        self.speclist = speclist

        self.stacked_spec = self.stacking()
        self.mf = None
        
        self.__init_check()
        

    def __init_check(self):
        if len(self.restfreqs) != len(self.speclist):
            functions.prt_info("Length of restfreqs and speclist should be equal.")
            return

    def __attr_check(self, primary=False, secondary=False):
        if len(self.speclist) == 0:
            functions.prt_info("No data is initialized for lines.speclist.")
            return
        
        for spec in self.speclist:
            if spec.__class__.__name__ != "obsspec":
                functions.prt_info("Invailed type in the lines.speclist, should be <class 'obsspec'>.")
                return 
        
        #only check the lines.speclist attribute
        if primary:
            return

        if self.stacked_spec is None:
            functions.prt_info("Line stacking is not performed.")
            functions.prt_info("Try: self.stacking(by='natural | rms')")
            return
        
        #check the lines.speclist and lines.stacked_spec attribute
        if secondary:
            return
            
        if self.mf is None:
            functions.prt_info("Matched filter is not performed.")
            functions.prt_info("Try: self.matched_filter(model)")
            return
    
    def stacking(self, by="natural"):
        '''
        line stacking for all the transitions.

        input: by="natural | rms"
        by: stacking methods.

        return: stacked_spec
        stacked_spec: 2d array of the stacked spectral
                      column one: velocity
                      column two: intensity of Ta
        '''
        self.__attr_check(primary=True)

        if len(self.restfreqs) == 0:
            return None
        
        if by not in ["natural", "rms"]:
            functions.prt_info("Unrecognized weighting method.")
            functions.prt_info("Only accept: 'natural' or 'weight'.")
            return None

        if by == "natural":
            N = len(self.restfreqs)
            velo_arr, int_arr = 0, 0
        
            for spec in self.speclist:
                velo_arr += spec.intensity[:, 1]
                int_arr += spec.intensity[:, 2]

            velo_arr /= N
            int_arr /= N

        if by == "rms":
            velo_arr, int_arr = [], []
            speclist_copy = np.copy(self.speclist)
            
            #find the peak intensity of each line
            peak_int_list = []
            for spec in speclist_copy:
                peak_int = np.nanmax(spec.intensity[:, 2])
                peak_int_list.append(peak_int)
            
            #maximum intensity of all lines
            max_int = np.nanmax(peak_int_list)

            #calculate the rms of the intensity for each line
            int_rms_list = []
            for eachline in range(len(self.restfreqs)):
                inten = speclist_copy[eachline].intensity[:, 2]

                rms = functions.get_rms(inten)

                #weighting by peak intensity and their rms of each line
                weight = peak_int_list[eachline] / max_int
                weight /= rms**2
                
                int_rms_list.append(rms)

                inten *= weight

                speclist_copy[eachline].intensity[:, 2] = inten

            #calculate the sum of the rms of each velocity step for all lines
            velo_rms_list = []
            velo_length = len(speclist_copy[0].intensity[:, 0])
            for velo in range(velo_length):
                sum_rms = 0
                for eachline in range(len(self.restfreqs)):
                    if np.isnan(speclist_copy[eachline].intensity[velo, 2]):
                        continue
                    else:
                        sum_rms += int_rms_list[eachline]**2
                velo_rms_list.append(sum_rms)

            velo_rms_list = np.array(velo_rms_list)
            velo_rms_list[velo_rms_list == 0] = np.nan

            for spec in speclist_copy:
                velo_arr.append(spec.intensity[:, 1])
                int_arr.append(spec.intensity[:, 2])

            velo_arr = np.array(velo_arr)
            int_arr = np.array(int_arr)
            
            velo_arr = np.mean(velo_arr, axis=0)
            int_arr = np.nansum(int_arr, axis=0) / velo_rms_list

        #the final stacked spectral
        stacked_spec = np.column_stack((velo_arr, int_arr))

        return stacked_spec

    def matched_filter(self, model):
        '''
        matched filter for the stacked line.

        input: model
        model: list of the gaussian line model.
               example: [peak_intensity, central_velocity, line_width]
        
        return: mf
        mf: 2d array of the matched filter.
            column one: velocity
            column two: SNR
        '''
        self.__attr_check(secondry=True)

        if type(model) != list or type(model) != type(np.array([])):
            functions.prt_info("Invailed input, must be list or array.")
            return None
        if len(model) != 3:
            functions.prt_info("Length of the array should be 3.")
            return None
        for i in model:
            if type(i) != float:
                functions.prt_info("Invailed input, must be float.")
                return None
            
        velo = self.stacked_spec[:, 0]
        amp, mu, sigma = model

        line_model = amp * np.exp(-0.5*((velo-mu)/sigma)**2)
        ss_int = self.stacked_spec[:, 1]

        mf = np.correlate(ss_int, line_model, mode="same")

        return np.column_stack((velo, mf))

    def show_spec(self, vpos=None, xunit="velo", overlay=False, tight=False, title=None):
        '''
        plot the intensity of all lines.

        input: overlay=False, tight=False

        return: none
        '''
        self.__attr_check(primary=True)

        #plot grid figures
        x_arr_list, y_arr_list, label_list = [], [], []
        
        for spec in self.speclist:
            if xunit == "freq":
                x_arr_list.append(spec.intensity[:, 0])
                vpos = self.restfreqs
            if xunit == "velo":
                x_arr_list.append(spec.intensity[:, 1])
            y_arr_list.append(spec.intensity[:, 2])
        
        for line in self.restfreqs:
            label_list.append(str(line) + "MHz")

        plotting.plot_spec(x_arr_list, y_arr_list, label_list, vpos=vpos, xunit=xunit, overlay=overlay, tight=tight, title=title)

    def show_ss(self, title=None):
        '''
        plot the stacked spectral.
        '''
        self.__attr_check(secondary=True)

        x_arr = self.stacked_spec[:, 0]
        y_arr = self.stacked_spec[:, 1]
        label = "stacked spectral"

        plotting.plot_spec([x_arr], [y_arr], [label], title=title)

    def show_mf(self, title=None):
        '''
        plot the matched filter.
        '''
        self.__attr_check()

        x_arr = self.mf[:, 0]
        y_arr = self.mf[:, 1]
        label = "SNR"

        plt.figure(figsize=(15,6))

        plt.plot(x_arr, y_arr, label=label)
        
        plt.title(title)
        plt.xlabel(r"V_lsr / km/s")
        plt.ylabel("SNR")
        plt.grid()
        plt.legend(loc="upper right")

        plt.show()

class molspec():
    '''
    molspec class to store all of the data.

    attributes: obssetup, name="", beam="M01", restfreqs=[], mjds=[], obsdates=[], lines=None, loadpath=None.
    obssetup: dict. some user defined key-value pairs.
    name: string. molecular name.
    beam: string. FAST beam.
    restfreqs: list of rest frequencies.
    mjds: list of mjd.
    obsdates: list of observational date.
    lines: lines class
    loadpath: string. read from files if provided.
              filename example: loadpath/obsdate/M01_Ta_1419_1421_ave.npy
    
    methods:
    load_from_file()
    '''
    def __init__(self, obssetup, name="", beam="M01", restfreqs=[], mjds=[], obsdates=[], lines=None, loadpath=None):
        self.obssetup = obssetup
        self.name = name
        self.beam = beam
        self.restfreqs = restfreqs
        self.mjds = mjds
        self.obsdates = obsdates
        self.lines = lines
        self.loadpath = loadpath
        
        self.__init_check()
        
        if self.loadpath != None:
            self.load_from_file()

    def __init_check(self):
        if type(self.obssetup) != dict:
            functions.prt_info("First argument should be type of dict.")
            return
        if "target_ra" not in self.obssetup:
            functions.prt_info("key: target_ra not in the dict.")
            return
        if "target_dec" not in self.obssetup:
            functions.prt_info("key: target_dec not in the dict.")
            return
        if "bandw" not in self.obssetup:
            functions.prt_info("key: bandw not in the dict.")
            return
        if "nchan" not in self.obssetup:
            functions.prt_info("key: nchan not in the dict.")
            return

    def load_from_file(self):
        if self.loadpath is None:
            functions.prt_info("loadpath is None. Nothing to do.")
            return
        if not functions.is_dir_exists(self.loadpath):
            functions.prt_info("path %s is not exist.", self.loadpath)
            return

        #init lines class
        self.lines = lines(restfreqs=[], speclist=[])
        self.lines.restfreqs = self.restfreqs
        
        #init self.lines.speclist
        nlines = len(self.restfreqs)
        for iline,eachline in enumerate(self.restfreqs):
            functions.prt_info("Processing for transition %f MHz for %s (%d/%d)", eachline, self.name, iline+1, nlines)
            
            #init obsspec class
            spec = obsspec(obsdates=[], mjds=[], specs=[])
            spec.obsdates = self.obsdates
            spec.mjds = self.mjds
            
            #init obsspec.specs. for each mjd spectral
            for eachmjd,eachobs in zip(self.mjds, self.obsdates):
                #choose the spectral data
                for fi in sorted(os.listdir(self.loadpath+"/"+eachobs)):
                    if (self.beam not in fi) or ("Ta" not in fi) or ("ave.npy" not in fi):
                        continue
                    
                    freqmin, freqmax = fi.split("_")[2:4]
                    freqmin, freqmax = int(freqmin), int(freqmax)

                    if eachline < freqmin or eachline > freqmax:
                        continue
                    
                    functions.prt_info("Initialization from data file %s...", eachobs+"/"+fi)
                    
                    #original frequency array
                    freq_arr = functions.get_freq_range(xmin=freqmin, xmax=freqmax)

                    #velocity correction
                    #vlsr = pipeline.vlsr_corr(self.obssetup["target_ra"], self.obssetup["target_dec"], eachmjd)

                    #doppler correction
                    #doppler = pipeline.rest_corr(vlsr, eachline)
                    
                    #print(eachline+doppler)

                    #correction of the original frequency array
                    #freq_arr = pipeline.freq_shift(freq_arr, doppler, resolution=self.obssetup["bandw"]/self.obssetup["nchan"])

                    #convertion from frequency domain to velocity domain
                    velo_arr = pipeline.freq_to_velo(freq_arr, eachline)

                    #shifting Ta
                    ta_arr = np.load(self.loadpath+"/"+eachobs+"/"+fi)
                    #ta_arr = pipeline.ta_shift(ta_arr, doppler, resolution=self.obssetup["bandw"]/self.obssetup["nchan"])
                    
                    #remove baseline
                    window = 0.3    #MHz
                    baseline = pipeline.medfilter(ta_arr, window)
                    
                    #xx = [freq_arr, freq_arr]
                    #yy = [ta_arr, baseline]
                    #ll = ["ta", "baseline"]
                    #plotting.plot_spec(xx, yy, ll, xunit="freq", overlay=True, title=eachobs)
                    #plotting.plot_spec([freq_arr], [ta_arr-baseline], ["baseline removed"], xunit="freq")
                    
                    ta_arr = ta_arr - baseline
                    
                    #generate the 2d array of the spectral data
                    xyz = np.column_stack((freq_arr, velo_arr, ta_arr))
                    
                    spec.specs.append(xyz)
                #end for each file
            #end for each mjd
            spec.intensity = spec.ave_spec()
            #spec.show_int()
            #spec.intensity = spec.int_cutoff(vmin=-40, vmax=40)
            self.lines.speclist.append(spec)
        #end for each line

    
    def xxx(self):
        pass
    

