#!/usr/bin/python3
#-*- encoding: utf-8 -*-
import sys
if sys.version_info.major < 3:
    print("This Code runs under Python 3.")
    exit(0)
import os
from itertools import product
from uwb_tool import pipeline
from uwb_tool import functions


def cal_comet_raw_data(comet_name, obsdate, receiver, freq_limit):
    obssetup = {}
    obssetup["comet"] = comet_name
    obssetup["keep_freq"] = freq_limit
    obssetup["receiver"] = receiver
    obssetup["bandw"] = functions.UWB_receiver[receiver]["band_width"]

    obssetup["sampling_time"] = 1
    obssetup["period_on"] = 1
    obssetup["period_off"] = 9
    obssetup["source_on"] = 300
    obssetup["source_switch"] = 30

    cal_file = functions.UWB_receiver[receiver]["cal_file"]
    obssetup["tcal_file"] = f"./tcal/20220617-high-cal/{cal_file}"
    obssetup["fits_files"] = f"comet_A3_user-defined_spec-{receiver}_*.fits"


    #### FAST cluster setup
    my_disk = "/media/longfei/c57a45a4-0626-4f7a-bdbb-f0d65a153c9d"
    base_input_path = f"{my_disk}/N2024_8/comet_A3/{obsdate}"
    base_output_path = f"{my_disk}/N2024_8_only_for_OH"
    src_path = base_input_path
    out_path = "/".join([base_output_path, comet_name, obsdate, receiver])

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


    pipeline.comet_cal_power_data(obssetup, power_source_on_path, tacal_source_on_path, freqlimit=freq_limit, cycle=6, fig_path=fig_source_on_path)
    pipeline.comet_cal_power_data(obssetup, power_source_off_path, tacal_source_off_path, freqlimit=freq_limit, cycle=6, fig_path=fig_source_off_path)

    pipeline.ta_ave(obssetup, tacal_source_on_path, prod_source_on_path, freqlimit=freq_limit)
    pipeline.ta_ave(obssetup, tacal_source_off_path, prod_source_off_path, freqlimit=freq_limit)
    pipeline.comet_ON_minus_OFF(obssetup, prod_path, freqlimit=freq_limit)

    pipeline.merge_figs(obssetup, fig_source_on_path, prod_path+"/sourceONfigs", freqlimit=freq_limit)
    pipeline.merge_figs(obssetup, fig_source_off_path, prod_path+"/sourceOFFfigs", freqlimit=freq_limit)


if __name__ == "__main__":
    comet_name = "A3"
    obsdate = ["20241003", "20241004", "20241013", "20241014", "20241016", "20241017"]
    
    receiver = "UWB2"
    freq_limit_list = [[1663,1669]]
    for obs_date,freq_limit in product(obsdate, freq_limit_list):
        cal_comet_raw_data(comet_name, obs_date, receiver, freq_limit)
    
    
    