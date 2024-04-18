#!/usr/bin/python3
#-*- encoding: utf-8 -*-
import sys
if sys.version_info.major < 3:
    print("This Code runs under Python 3.")
    exit(0)


import os
import numpy as np
from uwb_tool import pipeline
from uwb_tool import functions


obssetup = {}
obssetup["ra"] = None
obssetup["dec"] = None
obssetup["receiver"] = "UWB2"
obssetup["bandw"] = functions.UWB_receiver[obssetup["receiver"]]["band_width"]
obssetup["nchan"] = 1048576

obssetup["sampling_time"] = 1
obssetup["period_on"] = 1
obssetup["period_off"] = 9
obssetup["source_on"] = 20
obssetup["source_switch"] = 30

obssetup["tcal_file"] = "./tcal/20220617-high-cal/0.95-1.75GHz-Th.xlsx"
obssetup["fits_files"] = "Comet_UWB_user-defined-UWB2_*.fits"


# comet_name = "arp220"
# base_input_path = "./comet_uwb_obs_data"
# base_output_path = "./comet_outputs"
# src_path = "/".join([base_input_path, comet_name])
# out_path = "/".join([base_output_path, comet_name])

comet_name = "12P"
base_input_path = "/data31/N2023_9/Comet_UWB/20240417"
base_output_path = "./comet_outputs"
src_path = base_input_path
out_path = "/".join([base_output_path, comet_name])

pol_averaged_path = "/".join([out_path, "pol_averaged"])

source_onoff_split_path = "/".join([out_path, "source_onoff"])
power_source_on_path = source_onoff_split_path+"/"+"sourceON"
power_source_off_path = source_onoff_split_path+"/"+"sourceOFF"

tacal_data_path = "/".join([out_path, "tacal"])
tacal_source_on_path = tacal_data_path+"/"+"sourceON"
tacal_source_off_path = tacal_data_path+"/"+"sourceOFF"

fig_path = "/".join([out_path, "plots"])
fig_source_on_path = fig_path + "/" + "sourceON"
fig_source_off_path = fig_path + "/" + "sourceOFF"

prod_path = "/".join([out_path, "product"])
prod_source_on_path = prod_path + "/" + "sourceON"
prod_source_off_path = prod_path + "/" + "sourceOFF"


if not os.path.exists(prod_source_off_path + "/power_ave.npy"):
    pipeline.polave_fits(obssetup, src_path, pol_averaged_path, fig_path=fig_path)

    load_path = pol_averaged_path
    pipeline.comet_split_source_onoff(obssetup, load_path, source_onoff_split_path)

    pipeline.power_ave(power_source_on_path, prod_source_on_path)
    pipeline.power_ave(power_source_off_path, prod_source_off_path)


for freq_limit in [[1415, 1425], [1660, 1670]]:
    pipeline.comet_cal_power_data(obssetup, power_source_on_path, tacal_source_on_path, freqlimit=freq_limit, cycle=10, fig_path=fig_source_on_path)
    pipeline.comet_cal_power_data(obssetup, power_source_off_path, tacal_source_off_path, freqlimit=freq_limit, cycle=10, fig_path=fig_source_off_path)

    pipeline.ta_ave(obssetup, tacal_source_on_path, prod_source_on_path, freqlimit=freq_limit)
    pipeline.ta_ave(obssetup, tacal_source_off_path, prod_source_off_path, freqlimit=freq_limit)
    pipeline.comet_ON_minus_OFF(obssetup, prod_path, freqlimit=freq_limit)

    pipeline.merge_figs(obssetup, fig_source_on_path, prod_path+"/sourceONfigs", freqlimit=freq_limit)
    pipeline.merge_figs(obssetup, fig_source_off_path, prod_path+"/sourceOFFfigs", freqlimit=freq_limit)
