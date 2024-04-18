#!/usr/bin/python3
#-*- encoding: utf-8 -*-

'''
some plots
'''
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


from uwb_tool import functions

def plot_xy(xarr, yarr, xlabel="", ylabel="", title=""):
    plt.figure(figsize=(16,5))
    plt.grid()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    

    plt.plot(xarr, yarr)
    
    plt.show()

def plot_pol(fig_path, idx, xx, yy):
    '''
    input: fig_path, idx, xx, yy
    return: none

    saved file: fig_path + "/" + "_".join(["pol", "{}.pdf".format("%0.4d"%idx)])
    '''
    plt.figure(figsize=(16,5))
    plt.grid()
    plt.xlabel("Channel", fontsize=15)
    plt.ylabel("Power", fontsize=15)
    plt.title("polarization", fontsize=18)
    plt.yscale("log")

    plt.plot(xx, label="xx")
    plt.plot(yy, label="yy")

    plt.legend(loc="upper right", fontsize=15)

    fig_name = "_".join(["pol", "{}.pdf".format("%0.4d"%idx)])
    fig_name = fig_path + "/" + fig_name
    plt.savefig(fig_name)
    
    plt.close()
    del xx, yy

def plot_dyn(fig_path, freqlimit, idx, xyrange, image):
    '''
    input: fig_path, freqlimit, idx, xyrange, image
    return: none
    
    saved file: fig_path + "/" + "_".join([str(freqlimit[0]), str(freqlimit[1]), "%0.4d"%idx, "dyn.pdf"])
    '''
    plt.figure(figsize=(15,10))
    plt.xlabel("Frequency (MHz)", fontsize=15)
    plt.ylabel("Time (second)", fontsize=15)
    plt.xticks(size=13)
    plt.yticks(size=13)

    plt.imshow(image, aspect="auto", origin="lower", extent=xyrange)

    cb = plt.colorbar()
    cb.set_label("Power", fontsize=15)
    cb.ax.tick_params(labelsize=13)


    fig_name = "_".join([str(freqlimit[0]), str(freqlimit[1]), "%0.4d"%idx, "dyn.pdf"])
    fig_name = fig_path + "/" + fig_name
    plt.savefig(fig_name)
    
    plt.close()
    del xyrange, image

def plot_timebin(fig_path, freqlimit, idx, period, x, y):
    '''
    input: fig_path, freqlimit, idx, period, x, y
    return: none
    
    saved file: fig_path + "/" + "_".join([str(freqlimit[0]), str(freqlimit[1]), "%0.4d"%idx, "timebin.pdf"])
    '''
    plt.figure(figsize=(15,6))
    plt.grid()
    plt.xlabel("Time (second)", fontsize=15)
    plt.ylabel("Power", fontsize=15)
    plt.xticks(size=13)
    plt.yticks(size=13)

    plt.xlim((idx-1)*period, idx*period)

    plt.scatter(x, y)

    fig_name = "_".join([str(freqlimit[0]), str(freqlimit[1]), "%0.4d"%idx, "timebin.pdf"])
    fig_name = fig_path + "/" + fig_name
    plt.savefig(fig_name)
    
    plt.close()
    del x, y

def plot_Ta(fig_path, freqlimit, idx, freq_array, Ta_array):
    '''
    input: fig_path, freqlimit, idx, freq, Ta
    return: none
    
    saved file: fig_path + "/" + "_".join([str(freqlimit[0]), str(freqlimit[1]), "%0.4d"%idx, "Ta.pdf"])
    '''
    plt.figure(figsize=(15,6))
    plt.grid()
    plt.xlabel("Frequency (MHz)", fontsize=15)
    plt.ylabel("Ta (K)", fontsize=15)
    plt.xticks(size=13)
    plt.yticks(size=13)

    plt.xlim(freqlimit[0], freqlimit[1])

    plt.plot(freq_array, Ta_array)

    fig_name = "_".join([str(freqlimit[0]), str(freqlimit[1]), "%0.4d"%idx, "Ta.pdf"])
    fig_name = fig_path + "/" + fig_name
    plt.savefig(fig_name)
    
    plt.close()
    del freq_array, Ta_array


def plot_Cal(fig_path, freqlimit, idx, freq_array, cal_array, label):
    '''
    input: fig_path, freqlimit, idx, freq_array, cal_array, label
    return: none
    
    saved file: fig_path + "/" + "_".join([str(freqlimit[0]), str(freqlimit[1]), "%0.4d"%idx, "cal.pdf"])
    '''
    plt.figure(figsize=(15,6))
    plt.grid()
    plt.xlabel("Frequency (MHz)", fontsize=15)
    plt.ylabel("power", fontsize=15)
    plt.xticks(size=13)
    plt.yticks(size=13)

    plt.xlim(freqlimit[0], freqlimit[1])

    plt.plot(freq_array, cal_array)

    fig_name = "_".join([str(freqlimit[0]), str(freqlimit[1]), "%0.4d"%idx, f"{label}.pdf"])
    fig_name = fig_path + "/" + fig_name
    plt.savefig(fig_name)
    
    plt.close()
    del freq_array, cal_array

  
def plot_power(power_files, labels):
    '''
    plot the fullband data.

    input: power_files, labels
    power_files: list of string. file name
    labels: list of string

    return: none

    file example: M01_power_ave.npy
    interactive plot.
    '''
    power_files = list(power_files)
    labels = list(labels)

    if len(power_files) != len(labels):
        print("inputs should have same length.")
        return

    plt.figure(figsize=(15,6))
    plt.grid()
    plt.xlabel("Frequency (MHz)", fontsize=15)
    plt.ylabel("Power")
    plt.xticks(size=13)
    plt.yticks(size=13)

    freq = functions.get_freq_range(xmin=1000, xmax=1500, nchan=1048576, bandw=500)

    for fi,label in zip(power_files, labels):
        power = np.load(fi)
        plt.plot(freq, power, label=label)

    plt.legend(loc="upper right", fontsize=18)
    plt.show()
        

def plot_spec(x_arr_list, y_arr_list, label_list, vpos=None, xunit="velo", overlay=False, tight=False, title=None):
    '''
    grid plot of obsspec.specs (for different mjds) or lines.speclist.intensity (for different transitions) date.

    input: x_arr_list, y_arr_list, label_list, xunit="velo", overlay=False, tight=False
    x_arr_list: list of x-axis data
    y_arr_list: list of y-axis data
    label_list: list of labels
    xunit: type of the x-axis of the spectral. velocity or frequency
          example: xunit = "velo" or "freq"
    overlay: boolean
             True: stacking plot
             False: graid plot
    tight: boolean
           True: no gap between grid

    return: none
    '''
    subfigs = len(x_arr_list)
    
    col = 4 if subfigs > 4 else subfigs
    row = subfigs // col
    row += 0 if subfigs % col == 0 else 1
      
    width, height = 5, 3
    fig = plt.figure(figsize=(width*col,height*row))

    if overlay:
        ax = fig.add_subplot(1, 1, 1)
    elif tight:
        fig.subplots_adjust(wspace=0, hspace=0)
    # fig.subplots_adjust(left=0.074, right=0.987)

    for ii in range(row):
        for jj in range(col):
            cur_ax = col * ii + jj
            if cur_ax > subfigs - 1:
                break
            
            label = label_list[cur_ax]
            x = x_arr_list[cur_ax]
            y = y_arr_list[cur_ax]
            
            if not overlay:
                ax = fig.add_subplot(row, col, cur_ax+1)
            
            if vpos != None:
                if type(vpos) == float:
                    ax.axvline(vpos, ls="--", color="red")
                if type(vpos) == list:
                    ax.axvline(vpos[cur_ax], ls="--", color="red")
                
            if (vpos != None) and (xunit == "freq"):
                lr_width = 0.8
                ud_width = 0.08
                if type(vpos) == float:
                    ax.set_xlim(vpos-lr_width, vpos+lr_width)
                if type(vpos) == list:
                    ax.set_xlim(vpos[cur_ax]-lr_width, vpos[cur_ax]+lr_width)
                ax.set_ylim(-ud_width, +ud_width)
                
                #ax.set_xlim(vpos-0.3, vpos+0.8)
                #ax.set_ylim(-0.5, +9)


            ax.plot(x, y, label=label)

            ax.xaxis.grid()
            ax.yaxis.grid()
            ax.tick_params(axis='x', labelsize=13)
            ax.tick_params(axis='y', labelsize=13)

            # ax.ticklabel_format(style="plain", useOffset=False)
            # ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
            # ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
            # ax.yaxis.set_major_locator(ticker.MultipleLocator(5.0))
            # ax.yaxis.set_minor_locator(ticker.MultipleLocator(1.0))
            
            ax.legend(loc="upper left", fontsize=14)


            if not tight:
                if xunit == "velo":
                    ax.set_xlabel(r"V_${lsr}$ (km/s)", fontsize=13)
                if xunit == "freq":
                    ax.set_xlabel("Frequency (MHz)", fontsize=13)
                ax.set_ylabel("Ta (K)", fontsize=13)
                continue
                
            if ii == row-1 and jj == 0:
                if xunit == "velo":
                    ax.set_xlabel(r"V_${lsr}$ (km/s)", fontsize=13)
                if xunit == "freq":
                    ax.set_xlabel("Frequency (MHz)", fontsize=13)
                ax.set_ylabel("Ta (K)", fontsize=13)
            else:
                ax.xaxis.set_major_formatter(plt.NullFormatter())
                ax.yaxis.set_major_formatter(plt.NullFormatter())
        
        if cur_ax > subfigs - 1:
            break
    
    fig.suptitle(title, fontsize=20)
    plt.show()

