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
    obssetup["receiver"] = receiver
    obssetup["bandw"] = functions.UWB_receiver[receiver]["band_width"]

    obssetup["sampling_time"] = 1
    obssetup["period_on"] = 1
    obssetup["period_off"] = 9
    obssetup["source_on"] = 300
    obssetup["source_switch"] = 30

    cal_file = functions.UWB_receiver[receiver]["cal_file"]
    obssetup["tcal_file"] = f"./tcal/20220617-high-cal/{cal_file}"
    obssetup["fits_files"] = f"Comet_UWB_user-defined-{receiver}_*.fits"


    #### FAST cluster setup
    # base_input_path = f"/data31/N2023_9/Comet_UWB/{obsdate}"
    my_disk = "/media/longfei/c57a45a4-0626-4f7a-bdbb-f0d65a153c9d"
    base_input_path = f"{my_disk}/N2023_9/Comet_UWB/{obsdate}"
    base_output_path = f"{my_disk}/N2023_9_outputs"
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
    comet_name = "12P"
    obsdate = ["20240417", "20240424", "20240429", "20240503", "20240510", "20240511", "20240513"]
    
    receiver = "UWB1"
    freq_limit_list = [[500,650], [650,800], [800,950]]
    for obs_date,freq_limit in product(obsdate, freq_limit_list):
        cal_comet_raw_data(comet_name, obs_date, receiver, freq_limit)
    
    receiver = "UWB2"
    freq_limit_list = [[950,1150], [1150,1350], [1350,1550], [1550,1750]]
    for obs_date,freq_limit in product(obsdate, freq_limit_list):
        cal_comet_raw_data(comet_name, obs_date, receiver, freq_limit)
    
    receiver = "UWB3"
    freq_limit_list = [[1750,1950], [1950,2150], [2150,2350], [2350,2550]]
    for obs_date,freq_limit in product(obsdate, freq_limit_list):
        cal_comet_raw_data(comet_name, obs_date, receiver, freq_limit)
    
    receiver = "UWB4"
    freq_limit_list = [[2550,2700], [2700,2850], [2850,3000], [3000,3150], [3150,3300], [3300,3450]]
    for obs_date,freq_limit in product(obsdate, freq_limit_list):
        cal_comet_raw_data(comet_name, obs_date, receiver, freq_limit)
    
    