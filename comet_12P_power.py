#!/usr/bin/python3
#-*- encoding: utf-8 -*-
import sys
if sys.version_info.major < 3:
    print("This Code runs under Python 3.")
    exit(0)
import numpy as np
import matplotlib.pyplot as plt
from uwb_tool import functions

obssetup = {}
receiver = "UWB2"
obssetup["receiver"] = receiver
obssetup["bandw"] = functions.UWB_receiver[obssetup["receiver"]]["band_width"]


comet_name = "12P"
obsdate = "20240510"

power_file_path = f"./comet_outputs/{comet_name}/{obsdate}/product/sourceON/power_ave.npy"

power = np.load(power_file_path.format(comet_name, obsdate))

_, freq_array = functions.get_freq_mask_range(receiver)

plt.plot(freq_array, power)

plt.xlabel("Freqency (MHz)")
plt.ylabel("Ta (K)")
plt.grid()
plt.show()

