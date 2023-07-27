"""
Functions for XRD simulations
"""

import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
rc("text", usetex=True)
rc("font", **{"family":"sans-serif","sans-serif":["Helvetica"]},size="14")
rc("text.latex",preamble=r"\usepackage{sfmath}")

import os
from colour import Color
import Dans_Diffraction as df
import numpy as np

#------------------------------------------------------------------------------
''' simulates X-ray diffraction patterns '''

def full_XRD_sim(cif_dir, cif, wavelength, tt_max, pw, bg):
    #--------------------------------------------------------------------------
    # cif_dir (str) - file path to load cif from
    # cif (str) -- name of cif file
    # wavelength (float) -- instrument wavelength
    # tt_max (float) -- maximum 2theta value
    # pw (float) -- peak width, units of A^-1
    # bg (float) -- average of normal background
    
    # returns lists of calcualted Q and intensity values
    #--------------------------------------------------------------------------
    
    # load cif file as structure 
    path = os.path.join(cif_dir, cif + ".cif")
    struct = df.Crystal(path)
    
    # setup diffraction parameters
    energy_kev = df.fc.wave2energy(wavelength)
    struct.Scatter.setup_scatter("xray")
    wavevector_max = df.fc.calqmag(tt_max, energy_kev)
    
    # simultate XRD pattern
    q, ints = struct.Scatter.generate_powder(wavevector_max, 
                                             peak_width=pw, background=bg,
                                             powder_average=True)
    return q, ints


#-----------------------------------------------------------------------------
''' saves text file of simulated XRD pattern '''

def save_sim(directory, name, q, ints):
    #--------------------------------------------------------------------------
    # directory (str) -- file path to save data to
    # name (str) -- file name 
    # q (list) -- list of calculated Q values
    # ints (list) -- list of calculated intensity values
    #--------------------------------------------------------------------------
    
    # create file path
    path = os.path.join(directory, name + ".txt")
    
    # create new file
    new_file = open(path, "w")
    
    # write simulation data to new text file
    with new_file as f:
        for (q, ints) in zip(q, ints):
            f.write("{0} {1}\n".format(q, ints))
    new_file.close()
    
    
#-----------------------------------------------------------------------------
''' import XRD simulation text file '''

def import_sim(directory, name):
    #--------------------------------------------------------------------------
    # path (str) -- path to folder where text file is stored
    # name (str) -- text file name
    
    # returns lists of Q and intensity values from text file
    #--------------------------------------------------------------------------
    
    # create data file path 
    path = os.path.join(directory, name + ".txt")
    
    # load data from text file
    q, ints = np.loadtxt(path, unpack=True, dtype=float)
    
    return q, ints
    

#-----------------------------------------------------------------------------
''' load cif as structure '''

def load_cif(directory, name):
    #--------------------------------------------------------------------------
    # directory (str) -- path to folder where cif is stored
    # name (str) -- cif file name
    
    # returns structure as Crystal object
    #--------------------------------------------------------------------------
    
    # creates cif file path
    path = os.path.join(directory, name + ".cif")
    
    # convert cif to Crystal object
    struct = df.Crystal(path)
    
    return struct


#-----------------------------------------------------------------------------
''' set up diffraction simulation conditions '''

def sim_setup(struct, wl, max_tt):
    #--------------------------------------------------------------------------
    # struct -- Crystal structure object
    # wl (float) -- instrument wavelength
    # max_tt (float) -- maximum 2theta value
    
    # returns maximum wavevector for simulation
    #--------------------------------------------------------------------------
    
    energy_kev = df.fc.wave2energy(wl)
    struct.Scatter.setup_scatter("xray")
    max_wavevector = df.fc.calqmag(max_tt, energy_kev)
    
    return max_wavevector
    

#-----------------------------------------------------------------------------
''' simulate powder diffraction '''

def sim(struct, max_wv, pw, bg):
    #--------------------------------------------------------------------------
    # struct -- Crystal structure object
    # max_wv (float) -- maximum wavevector
    # pw (float) -- peak width, units of A^-1
    # bg (float) -- average of normal background
    
    # returns lists of calcualted Q and intensity values
    #--------------------------------------------------------------------------
    
    q, intensity = struct.Scatter.generate_powder(max_wv, peak_width=pw, 
                                                  background=bg, 
                                                  powder_average=True)
    return q, intensity
    

#-----------------------------------------------------------------------------
''' normalize intensity values '''

def norm(ints):
    #--------------------------------------------------------------------------
    # ints -- list of calculated intensity values
    
    # returns list of normalilzed intensity values
    #--------------------------------------------------------------------------
    
    norm_ints = (ints / np.max(ints))
    
    return norm_ints


#-----------------------------------------------------------------------------
''' generate (hkl) values '''

def hkl_sim(cif_dir, cif_name, wavelength, tt_max):
    #--------------------------------------------------------------------------
    # cif_dir (str) - file path to load cif from
    # cif_name (str) -- name of cif file
    # wavelength (float) -- instrument wavelength
    # tt_max (float) -- maximum 2theta value
    
    # returns (hkl), 2theta, and intensity as text
    #--------------------------------------------------------------------------
    
    # load cif file as structure
    path = os.path.join(cif_dir, cif_name + ".cif")
    struct = df.Crystal(path)
    
    # setup diffraction parameters
    e_kev = df.fc.wave2energy(wavelength)
    struct.Scatter.setup_scatter(scattering_type="xray", 
                                 energy_kev=e_kev, 
                                 min_twotheta=0, max_twotheta=tt_max)
    
    # generate (hkl) reflections
    reflections = struct.Scatter.print_all_reflections(min_intensity=0.1)
    
    return reflections


#------------------------------------------------------------------------------
''' save (hkl) values as text file '''

def save_hkl(directory, name, data):
    #--------------------------------------------------------------------------
    # directory (str) -- file path to save data to
    # name (str) -- file name 
    # data (str) -- (hkl) information
    #--------------------------------------------------------------------------
    
    # create file path
    path = os.path.join(directory, name + "_hkl.txt")
    
    # create new file
    new_file = open(path, "w")
    
    # write (hkl) data to new text file
    with new_file as f:
        f.write(data)
    new_file.close()
    
    
#------------------------------------------------------------------------------
''' calculate a difference curve '''
def diff_curve(q1, q2, ints1, ints2): 
    #--------------------------------------------------------------------------
    # q1 (nparray) -- x data of data set 1
    # q2 (nparray) -- x data of data set 2
    # ints1 (nparray) -- y data of data set 1
    # ints2 (nparray) -- y data of data set 2
    
    # returns array of x values and array of y values of difference curve
    #--------------------------------------------------------------------------
    
    diff_q_list = []
    diff_ints_list = []
    for i in range(0, len(q1)):
        q1_val = float('%.3f'%(q1[i]))
        for j in range(0, len(q2)):
            q2_val = float('%.3f'%(q2[j]))
            if q1_val == q2_val:
                diff_q_list.append(q1_val)
                diff_ints_list.append(ints1[i]-ints2[j])
    diff_q = np.array(diff_q_list)
    diff_ints = np.array(diff_ints_list)
    
    return diff_q, diff_ints


#------------------------------------------------------------------------------
''' save difference curve to text file '''
def save_diff_curve(q, ints, directory, name):
    #--------------------------------------------------------------------------
    # q (nparray) -- x data of difference curve
    # ints (nparray) -- y data of difference curve
    # directory (str) -- file path to save data to
    # name (str) -- file name
    #--------------------------------------------------------------------------
    
    # create file path
    path = os.path.join(directory, name + ".txt")
    
    # create new file
    new_file = open(path, "w")
    
    # write difference curve data to new text file
    with new_file as f:
        for (q, ints) in zip(q, ints):
            f.write("{0} {1}\n".format(q, ints))
    new_file.close()
    
    
#------------------------------------------------------------------------------
''' PLOTTING FUNCTIONS '''
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
''' generates a 2D color gradient '''

def gradient_gen_2D(upperleft, upperright, lowerleft, lowerright, rows, cols):
    #--------------------------------------------------------------------------
    # upperleft (str) -- hex code for upper left gradient color
    # upperright (str) -- hex code for upper right gradient color 
    # lowerleft (str) -- hex code for lower left gradient color 
    # lowerright (str) -- hex code for lower right gradient color
    # rows (int) -- number of rows
    # cols (int) -- number of columns
    
    # returns 2D array of hex codes
    #--------------------------------------------------------------------------
    
    # get colors from hex codes
    tl = Color(upperleft)
    tr = Color(upperright)
    bl = Color(lowerleft)
    br = Color(lowerright)
    
    # generate gradient from upper left to lower left
    start_col_colors = list(tl.range_to(bl, rows))
    
    # generate gradient from upper right to lower lright
    end_col_colors = list(tr.range_to(br, rows))
    
    # generate remaining gradient colors
    color_list = []
    for i in range(0, rows):
        c1 = start_col_colors[i]
        c2 = end_col_colors[i]
        curr_row_colors = list(c1.range_to(c2, cols))
        color_list.append(curr_row_colors)
    
    return color_list


#------------------------------------------------------------------------------
''' generate stacked plot of simulated data at fault probabilities 10, 20, 30, 
and 40%'''
def sim_4prob_stack(expt_q, expt_ints, p10_q, p10_ints, p20_q, p20_ints, p30_q, 
                   p30_ints, p40_q, p40_ints, y_min, y_max, title, fn=None, 
                   save_dir=None, saveFig=False):
    #--------------------------------------------------------------------------
    # expt_q (nparray) -- x data of experimental data set
    # expt_ints (nparray) -- y data of experimental data set
    # p10_q (nparray) -- x data of P=10% data set
    # p10_ints (nparray) -- y data of P=10% data set
    # p20_q (nparray) -- x data of P=20% data set
    # p20_ints (nparray) -- y data of P=20% data set
    # p30_q (nparray) -- x data of P=30% data set
    # p30_ints (nparray) -- y data of P=30% data set
    # p40_q (nparray) -- x data of P=40% data set
    # p40_ints (nparray) -- y data of P=40% data set
    # y_min (float) -- minimum y-axis value
    # y_max (float) -- maximum y-axis value
    # title (str) -- plot title
    # fn (str, optional) -- file name
    # save_dir (str, optional) -- file path to save plot to
    # saveFig (bool, optional) -- set to True to save .png to save_dir
    #--------------------------------------------------------------------------
    
    fig, (p) = plt.subplots(4, 1, figsize=(6,6))
    
    expt_min = np.min(expt_ints)
    
    q_data = [p10_q, p20_q, p30_q, p40_q]
    ints_data = [p10_ints, p20_ints, p30_ints, p40_ints]
    prob_list = [r"$P=10\%$", r"$P=20\%$", r"$P=30\%$", r"$P=40\%$"]
    colors = ["#00C6BF", "#009AE1", "#5D7AD3", "#B430C2"]
    
    # plot data
    for i in range(0, 4):
        p[i].scatter(expt_q, expt_ints, color="black", label="Observed", 
                     marker=".", s=8)
        p[i].plot(q_data[i], ints_data[i]+expt_min, color=colors[i], 
                  label="_" + prob_list[i], linewidth="2")
        
    # axis limits
    x_ax = (1.06, 1.85)
    y_ax = (y_min, y_max) 
        
    for i in range(0, 4):
        p[i].set_xlim(x_ax[0], x_ax[1])
        p[i].set_ylim(y_ax[0], y_ax[1])
        p[i].tick_params(axis="both", labelsize="12")
        if i <= 2:
            p[i].get_xaxis().set_visible(False)   
        
    # labels
    fig.supxlabel("Q (\AA" r"$^{-1}$, $\lambda=0.459744$" ")", fontsize="14")
    fig.supylabel("Intensity (counts, normalized)", fontsize="14")
    fig.suptitle(title, fontsize="14", y=1)
    
    peak_labels = [("(020)", 1.14), 
             ("(110)", 1.19), 
             ("(11-1)", 1.33), 
             ("(021)", 1.549), 
             ("(111)", 1.805)]

    for i in range(0, len(peak_labels)):
        p[0].text(peak_labels[i][1], y_ax[1]+0.005, peak_labels[i][0], 
                  ha="center", va="bottom", rotation=45, fontsize="10", 
                  color = "#4C4C4C")
    
    for i in range(0, 4):
        p[i].text(x_ax[1]-0.01, y_ax[1]-0.02, prob_list[i], fontsize="14", 
                  ha="right", va="top", color=colors[i])
        
    plt.subplots_adjust(hspace=0.1) 

    p[0].legend(handlelength=0.25, fontsize="12", loc="upper center")
    
    if saveFig == True:
        plt.savefig(save_dir + fn + "_sim_prob_stack" + ".png", 
                    bbox_inches="tight", pad_inches=0.5, dpi=1000) 


#------------------------------------------------------------------------------
''' generate stacked plot of differences of simulated and experimental data 
at fault probabilities 10, 20, 30, and 40%'''
def diff_4prob_stack(expt_q, expt_ints, p10_diff_q, p10_diff_ints, p20_diff_q, 
                     p20_diff_ints, p30_diff_q, p30_diff_ints, p40_diff_q, 
                     p40_diff_ints, y_min, y_max, title, fn=None, 
                     save_dir=None, saveFig=False):
    #--------------------------------------------------------------------------
    # expt_q (nparray) -- x data of experimental data set
    # expt_ints (nparray) -- y data of experimental data set
    # p10_diff_q (nparray) -- x data of expt vs P=10% difference data set
    # p10_diff_ints (nparray) -- y data of expt vs P=10% difference data set
    # p20_diff_q (nparray) -- x data of expt vs P=20% difference data set
    # p20_diff_ints (nparray) -- y data of expt vs P=20% difference data set
    # p30_diff_q (nparray) -- x data of expt vs P=30% difference data set
    # p30_diff_ints (nparray) -- y data of expt vs P=30% difference data set
    # p40_diff_q (nparray) -- x data of expt vs P=40% difference data set
    # p40_diff_ints (nparray) -- y data of expt vs P=40% difference data set
    # y_min (float) -- minimum y-axis value
    # y_max (float) -- maximum y-axis value
    # title (str) -- plot title
    # fn (str, optional) -- file name
    # save_dir (str, optional) -- file path to save plot to
    # saveFig (bool, optional) -- set to True to save .png to save_dir
    #--------------------------------------------------------------------------
    
    fig, (diff) = plt.subplots(4, 1, figsize=(6,6))
    
    expt_min = np.min(expt_ints)
    
    q_data = [p10_diff_q, p20_diff_q, p30_diff_q, p40_diff_q]
    ints_data = [p10_diff_ints, p20_diff_ints, p30_diff_ints, p40_diff_ints]
    prob_list = [r"$P=10\%$", r"$P=20\%$", r"$P=30\%$", r"$P=40\%$"]
    colors = ["#00C6BF", "#009AE1", "#5D7AD3", "#B430C2"]
    
    # plot
    for i in range(0, 4):
        diff[i].plot(q_data[i], ints_data[i]+expt_min, color="#BEBEBE", 
                     label=prob_list[i], linewidth="1.5")
        
    # axis limits
    x_ax = (1.06, 1.85)
    y_ax = (y_min, y_max) 
        
    for i in range(0, 4):
        diff[i].set_xlim(x_ax[0], x_ax[1])
        diff[i].set_ylim(y_ax[0], y_ax[1])
        diff[i].tick_params(axis="both", labelsize="12")
        if i <= 2:
            diff[i].get_xaxis().set_visible(False)   
        
    # labels
    fig.supxlabel("Q (\AA" r"$^{-1}$, $\lambda=0.459744$" ")", fontsize="14")
    ylabel = r"I$_{\mbox{\small Obs}} -$ I$_{\mbox{\small Calc}}$ (counts)"
    fig.supylabel(ylabel, fontsize="14")
    fig.suptitle(title, fontsize="14", y=1)
    
    peak_labels = [("(020)", 1.14), 
             ("(110)", 1.19), 
             ("(11-1)", 1.33), 
             ("(021)", 1.549), 
             ("(111)", 1.805)]

    for i in range(0, len(peak_labels)):
        diff[0].text(peak_labels[i][1], y_ax[1]+0.005, peak_labels[i][0], 
                     ha="center", va="bottom", rotation=45, fontsize="10", 
                     color = "#4C4C4C")
    
    for i in range(0, 4):
        diff[i].text(x_ax[1]-0.01, y_ax[1]-0.02, prob_list[i], fontsize="14", 
                     ha="right", va="top", color=colors[i])
        
    plt.subplots_adjust(hspace=0.1)
    
    if saveFig == True:
        plt.savefig(save_dir + fn + "_diff_prob_stack" + ".png", 
                    bbox_inches="tight", pad_inches=0.5, dpi=1000) 


#------------------------------------------------------------------------------
''' generate stacked simulated data and difference with experimental data at 
fault probabilities 10, 20, 30, and 40% '''
def sim_and_diff_prob_stack(expt_q, expt_ints, p10_q, p10_ints, p20_q, 
                            p20_ints, p30_q, p30_ints, p40_q, p40_ints, 
                            p10_diff_q, p10_diff_ints, p20_diff_q, 
                            p20_diff_ints, p30_diff_q, p30_diff_ints, 
                            p40_diff_q, p40_diff_ints, y_min, y_max, title, 
                            fn=None, save_dir=None, saveFig=False):
    #--------------------------------------------------------------------------
    # expt_q (nparray) -- x data of experimental data set
    # expt_ints (nparray) -- y data of experimental data set
    # p10_q (nparray) -- x data of P=10% data set
    # p10_ints (nparray) -- y data of P=10% data set
    # p20_q (nparray) -- x data of P=20% data set
    # p20_ints (nparray) -- y data of P=20% data set
    # p30_q (nparray) -- x data of P=30% data set
    # p30_ints (nparray) -- y data of P=30% data set
    # p40_q (nparray) -- x data of P=40% data set
    # p40_ints (nparray) -- y data of P=40% data set
    # p10_diff_q (nparray) -- x data of expt vs P=10% difference data set
    # p10_diff_ints (nparray) -- y data of expt vs P=10% difference data set
    # p20_diff_q (nparray) -- x data of expt vs P=20% difference data set
    # p20_diff_ints (nparray) -- y data of expt vs P=20% difference data set
    # p30_diff_q (nparray) -- x data of expt vs P=30% difference data set
    # p30_diff_ints (nparray) -- y data of expt vs P=30% difference data set
    # p40_diff_q (nparray) -- x data of expt vs P=40% difference data set
    # p40_diff_ints (nparray) -- y data of expt vs P=40% difference data set
    # y_min (float) -- minimum y-axis value
    # y_max (float) -- maximum y-axis value
    # title (str) -- plot title
    # fn (str, optional) -- file name
    # save_dir (str, optional) -- file path to save plot to
    # saveFig (bool, optional) -- set to True to save .png to save_dir
    #--------------------------------------------------------------------------

    fig, (p) = plt.subplots(2, 2, figsize=(7,7))
    
    expt_min = np.min(expt_ints)
    
    q_data = [p10_q, p20_q, p30_q, p40_q]
    ints_data = [p10_ints, p20_ints, p30_ints, p40_ints]
    diff_q_data = [p10_diff_q, p20_diff_q, p30_diff_q, p40_diff_q]
    diff_ints_data = [p10_diff_ints, p20_diff_ints, p30_diff_ints, 
                      p40_diff_ints]
    prob_list = [r"$P=10\%$", r"$P=20\%$", r"$P=30\%$", r"$P=40\%$"]
    colors = ["#00C6BF", "#009AE1", "#5D7AD3", "#B430C2"]
    
    # plot
    for i in range(0, 2):
        for j in range(0, 2):
            if i == 0:
                index = j
            elif i == 1:
                index = j+2
                
            p[i,j].scatter(expt_q, expt_ints, color="black", label="Observed", 
                           marker=".", s=8)
            p[i,j].plot(q_data[index], ints_data[index]+expt_min, 
                        color=colors[index], label="_" + prob_list[index], 
                        linewidth="2")
            p[i,j].plot(diff_q_data[index], diff_ints_data[index]-0.17, 
                        color="#BEBEBE", label="Difference", linewidth="1")
    
    # axis limits
    x_ax = (1.06, 1.85)
    y_ax = (y_min, y_max) 
        
    for i in range(0, 2):
            for j in range(0, 2):
                p[i,j].set_xlim(x_ax[0], x_ax[1])
                p[i,j].set_ylim(y_ax[0], y_ax[1])
                p[i,j].tick_params(axis="both", labelsize="12")
            if i == 0:
                p[i,j].get_xaxis().set_visible(False)
            if j == 1:
                p[i,j].get_yaxis().set_visible(False)
                
    ticks = p[1,0].xaxis.get_major_ticks()
    ticks[1].set_visible(False)
        
    # labels
    fig.supxlabel("Q (\AA" r"$^{-1}$, $\lambda=0.459744$" ")", fontsize="14", 
                  y=0.03)
    fig.supylabel("Intensity (counts, normalized)", fontsize="14")
    fig.suptitle(title, fontsize="14", y=0.93)
    
    for i in range(0, 2):
        for j in range(0,2):
            if i == 0:
                index = j
            elif i == 1:
                index = j+2
            
            p[i,j].text(x_ax[1]-0.01, y_ax[1]-0.02, prob_list[index], 
                        fontsize="14", ha="right", va="top", 
                        color=colors[index])
        
    plt.subplots_adjust(wspace=0.03, hspace=0.03) 

    p[1,1].legend(handlelength=1, fontsize="12", loc="lower right")
    
    if saveFig == True:
        plt.savefig(save_dir + fn + "_sim_and_diff_prob_stack" + ".png", 
                    bbox_inches="tight", pad_inches=0.5, dpi=1000) 


#------------------------------------------------------------------------------
''' generate plot of unfaulted supercell vs expt and faulted supercell vs 
expt'''
def ideal_v_flt(expt_q, expt_ints, id_q, id_ints, flt_q, flt_ints, id_diff_q, 
                id_diff_ints, flt_diff_q, flt_diff_ints, y_min, y_max, title, 
                prob_percent, s_vector_title, id_peak_y_pos, flt_peak_y_pos, 
                fn=None, save_dir=None, saveFig=False):
    #--------------------------------------------------------------------------
    # expt_q (nparray) -- x data of experimental data set
    # expt_ints (nparray) -- y data of experimental data set
    # id_q (nparray) -- x data of unfaulted supercell data set
    # id_ints (nparray) -- y data of unfaulted supercell data set
    # flt_q (nparray) -- x data of faulted supercell data set
    # flt_ints (nparray) -- y data of faulted supercell data set
    # id_diff_q (nparray) -- x data of expt vs unfaulted difference data set
    # id_diff_ints (nparray) -- y data of expt vs unfaulted difference data set
    # flt_diff_q (nparray) -- x data of expt vs faulted difference data set
    # flt_diff_ints (nparray) -- y data of expt vs faulted difference data set
    # y_min (float) -- minimum y-axis value
    # y_max (float) -- maximum y-axis value
    # title (str) -- plot title
    # prob_percent (int) -- fault probability percentage
    # s_vector_title (str) -- title information for faulted plot
    # id_peak_y_pos (list) -- list of y values for unfaulted plot peak labels
    # flt_peak_y_pos (list) -- list of y values for faulted plot peak labels
    # fn (str, optional) -- file name
    # save_dir (str, optional) -- file path to save plot to
    # saveFig (bool, optional) -- set to True to save .png to save_dir
    #--------------------------------------------------------------------------
    
    fig, (p) = plt.subplots(1, 2, sharey=True, figsize=(14,7))
    
    colors = ["#FA008E", "#2FF8B9"]
    
    # ideal plot
    p[0].scatter(expt_q, expt_ints, color="black", label="Observed", 
                 marker=".")
    p[0].plot(id_q, id_ints, color=colors[0], label="Calculated", 
              linewidth="2")
    p[0].plot(id_diff_q, id_diff_ints-0.2, color="#BEBEBE", label="Difference", 
              linewidth="1")

    # faulted plot
    p[1].scatter(expt_q, expt_ints, color="black", label="Observed", 
                 marker=".")
    p[1].plot(flt_q, flt_ints, color=colors[1], label="Calculated", 
              linewidth="2")
    p[1].plot(flt_diff_q, flt_diff_ints-0.2, color="#BEBEBE", 
              label="Difference", linewidth="1")
    
    # axis limits
    x_ax = (1.06, 1.85)
    y_ax = (y_min, y_max)
    
    for i in range(0, 2):
        p[i].set_xlim(x_ax[0], x_ax[1])
        p[i].set_ylim(y_ax[0], y_ax[1])
        
    p[1].get_yaxis().set_visible(False)

    # labels
    fig.supxlabel("Q (\AA" r"$^{-1}$, $\lambda=0.459744$" ")")
    fig.supylabel("Intensity (counts, normalized)", x=0.065)
    
    props = dict(boxstyle="round", facecolor="lavender", alpha=0.5, pad=0.25)
    id_title = "\n".join(("Unfaulted Supercell", r"$N = 1000$"))
    flt_title = "\n".join((str(prob_percent) + r"\% Fault Probability, " + 
                           s_vector_title, r"$N = 1000$"))
    
    p[0].text(x_ax[0]+0.02, y_ax[1]-0.025, id_title, bbox=props, ha="left", 
              va="top")
    p[1].text(x_ax[0]+0.02, y_ax[1]-0.02, flt_title, bbox=props, ha="left", 
              va="top")
    
    fig.suptitle(title, y=0.93)
    
    id_labels = [("(020)", 1.14, id_peak_y_pos[0]), 
                 ("(110)", 1.19, id_peak_y_pos[1]), 
                 ("(11-1)", 1.33, id_peak_y_pos[2]), 
                 ("(021)", 1.549, id_peak_y_pos[3]), 
                 ("(111)", 1.805, id_peak_y_pos[4])]
    
    flt_labels = [("(020)", 1.14, flt_peak_y_pos[0]), 
                 ("(110)", 1.19, flt_peak_y_pos[1]), 
                 ("(11-1)", 1.33, flt_peak_y_pos[2]), 
                 ("(021)", 1.549, flt_peak_y_pos[3]), 
                 ("(111)", 1.805, flt_peak_y_pos[4])]

    for i in range(0, len(id_labels)):
        p[0].text(id_labels[i][1], id_labels[i][2], id_labels[i][0], 
                  ha="center", va="center")
        p[1].text(flt_labels[i][1], flt_labels[i][2], flt_labels[i][0], 
                  ha="center", va="center")

    plt.subplots_adjust(wspace=0.05)

    p[0].legend(handlelength=1)
    p[1].legend(handlelength=1)

    if saveFig == True:
        plt.savefig(save_dir + fn + "_ideal_v_fault" + ".png", 
                    bbox_inches="tight", pad_inches=0.5, dpi=1000)


#------------------------------------------------------------------------------
''' generate plot of differences of unfaulted vs faulted simulated against
experimental data'''
def ideal_v_fault_fit_diff(expt_ints, id_diff_q, id_diff_ints, flt_diff_ints, 
                           y_min, y_max, fn=None, save_dir=None, 
                           saveFig=False):
    #--------------------------------------------------------------------------
    # expt_ints (nparray) -- y data of experimental data set
    # id_diff_q (nparray) -- x data of expt vs unfaulted difference data set
    # id_diff_ints (nparray) -- y data of expt vs unfaulted difference data set
    # flt_diff_ints (nparray) -- y data of expt vs faulted difference data set
    # y_min (float) -- minimum y-axis value
    # y_max (float) -- maximum y-axis value
    # fn (str, optional) -- file name
    # save_dir (str, optional) -- file path to save plot to
    # saveFig (bool, optional) -- set to True to save .png to save_dir
    #--------------------------------------------------------------------------

    plt.figure(figsize=(6, 1))
    
    expt_min = np.min(expt_ints)

    plt.plot(id_diff_q, (id_diff_ints - flt_diff_ints - expt_min), 
             color="#BEBEBE")

    # axis limits
    plt.xlim(1.06, 1.85)
    plt.xticks(fontsize=12)

    y_ax = (y_min, y_max)
    plt.ylim(y_ax[0], y_ax[1])
    plt.yticks(fontsize=12)

    # axis labels
    plt.xlabel("Q (\AA" r"$^{-1}$, $\lambda=0.459744$" ")", fontsize="14")
    ylabel = r"I$_{\mbox{\small UF}} -$ I$_{\mbox{\small F}}$" + "\n (counts)"
    plt.ylabel(ylabel, rotation="horizontal", ha="center", va="center", 
               labelpad=20)
    
    diff_labels = [("(020)", 1.14), 
                 ("(110)", 1.19), 
                 ("(11-1)", 1.33), 
                 ("(021)", 1.549), 
                 ("(111)", 1.805)]

    for i in range(0, len(diff_labels)):
        plt.text(diff_labels[i][1], y_ax[1]+0.005, diff_labels[i][0], 
                 ha="center", va="bottom", rotation=45, fontsize="12", 
                 color = "#BEBEBE")

    plt.subplots_adjust(hspace=0.5)    
    
    if saveFig == True:
        plt.savefig(save_dir + fn + "_ideal_v_fault_fit_diff" + ".png", 
                    bbox_inches="tight", pad_inches=0.5, dpi=1000)


#------------------------------------------------------------------------------
''' generate plot comparing fit improvement of different simulation models at 
fault probabilities 10, 20, 30, and 40% '''
def model_comparison(id_v_p10_q, id_v_p10_ints, id_v_p20_q, id_v_p20_ints, 
                     id_v_p30_q, id_v_p30_ints, id_v_p40_q, id_v_p40_ints,
                     y_min, y_max, fn=None, save_dir=None, saveFig=False):
    #--------------------------------------------------------------------------
    # id_v_p10_q (nparray) -- x data of unfaulted vs P=10% difference
    # id_v_p10_ints (nparray) -- y data of unfaulted vs P=10% difference
    # id_v_p20_q (nparray) -- x data of unfaulted vs P=20% difference
    # id_v_p20_ints (nparray) -- y data of unfaulted vs P=20% difference
    # id_v_p30_q (nparray) -- x data of unfaulted vs P=30% difference
    # id_v_p30_ints (nparray) -- y data of unfaulted vs P=30% difference
    # id_v_p40_q (nparray) -- x data of unfaulted vs P=40% difference
    # id_v_p40_ints (nparray) -- y data of unfaulted vs P=40% difference
    # y_min (float) -- minimum y-axis value
    # y_max (float) -- maximum y-axis value
    # fn (str, optional) -- file name
    # save_dir (str, optional) -- file path to save plot to
    # saveFig (bool, optional) -- set to True to save .png to save_dir
    #--------------------------------------------------------------------------
    
    fig, (p) = plt.subplots(4, 5, figsize=(9, 5))
    
    gradient = gradient_gen_2D("#00C6BF", "#009AE1", "#5D7AD3", "#B430C2", 4,5)
    prob_list = [r"$P=10\%$", r"$P=20\%$", r"$P=30\%$", r"$P=40\%$"]
    
    # (peak name, x_min, x_max)
    peaks = [("(020)", 1.1305, 1.1495),
            ("(110)", 1.1805, 1.1995),
            ("(11-1)", 1.3205, 1.3395),
            ("(021)", 1.5405, 1.5595),
            ("(111)", 1.7905, 1.8095)]
    
    for col in range(0, 5):
        p[0,col].plot(id_v_p10_q, id_v_p10_ints, color=gradient[0][col].hex)
        p[1,col].plot(id_v_p20_q, id_v_p20_ints, color=gradient[1][col].hex)
        p[2,col].plot(id_v_p30_q, id_v_p30_ints, color=gradient[2][col].hex)
        p[3,col].plot(id_v_p40_q, id_v_p40_ints, color=gradient[3][col].hex)
        
    for col in range(0, 5):
        for row in range(0, 4):
            p[row,col].set_xlim(peaks[col][1], peaks[col][2])
            p[row,col].set_ylim(y_min, y_max)
            if row != 3:
                p[row,col].get_xaxis().set_visible(False)
            if row == 3:
                p[row,col].tick_params(axis="x", labelsize="10", rotation=30)
            if col != 0:
                p[row,col].get_yaxis().set_visible(False)
            if col == 0:
                p[row,col].tick_params(axis="y", labelsize="10")
                
    fig.supxlabel("Q (\AA" r"$^{-1}$, $\lambda=0.459744$" ")", fontsize="14", 
                  y=-0.03)
    ylabel = r"I$_{\mbox{\small UF}} -$ I$_{\mbox{\small F}}$ (counts)"
    fig.supylabel(ylabel, fontsize="14", x=0.05)
    
    for col in range(0, 5):
        x_mid = ((peaks[col][2] - peaks[col][1]) / 2) + peaks[col][1]
        p[0,col].text(x_mid, y_max, peaks[col][0] + " Reflection", ha="center", 
                      va="bottom", fontsize="12", color=gradient[0][col].hex)
        
    for row in range(0, 4):
        y_mid = ((y_max - y_min) / 2) + y_min
        p[row,3].text(1.581, y_mid, prob_list[row], ha="left", va="center", 
                      fontsize="12", color=gradient[row][3].hex)
        
        plt.subplots_adjust(hspace=0.1, wspace=0.1)
    
    if saveFig == True:
        plt.savefig(save_dir + fn + "_model_compare" + ".png", 
                    bbox_inches="tight", pad_inches=0.5, dpi=1000)






