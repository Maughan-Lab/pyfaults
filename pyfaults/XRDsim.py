'''
Functions for XRD simulations
'''

from colour import Color
import Dans_Diffraction as df
import numpy as np
import re

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.pyplot import rc
rc("text", usetex=True)
rc("font", **{"family":"sans-serif","sans-serif":["Helvetica"]},size="14")
rc("text.latex",preamble=r"\usepackage{sfmath}")

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' XRD SIMULATION FUNCTIONS '''
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def fullSim(path, cif, wl, tt_max, pw=None, bg=None):
    '''
    Simulates XRD pattern

    Parameters
    ----------
    path : str
        File directory where cif is store
    cif : str
        File name
    wl : float
        Instrument wavelength (A)
    tt_max : float
        Maximum 2theta (deg)
    pw : float, optional
        Peak width (A^-1)
    bg : float, optional
        Average of normal background

    Returns
    -------
    q : array
        Calculated Q values
    ints : array
        Calculated intensity values

    '''
    
    # load cif
    struct = df.Crystal(path + cif + ".cif")
    
    # setup diffraction parameters
    energy_kev = df.fc.wave2energy(wl)
    struct.Scatter.setup_scatter("xray")
    wavevector_max = df.fc.calqmag(tt_max, energy_kev)
    
    if pw is None:
        pw = 0.0
    if bg is None:
        bg = 0
    
    # simultate XRD pattern
    q, ints = struct.Scatter.generate_powder(wavevector_max, peak_width=pw, 
                                             background=bg, powder_average=True)
    return q, ints


#-----------------------------------------------------------------------------
def saveSim(path, fn, q, ints):
    '''
    Saves simulated XRD as text file

    Parameters
    ----------
    path : str
        File directory
    fn : str
        File name
    q : array
        Calculated Q values
    ints : array
        Calculated intensity values

    Returns
    -------
    None

    '''
    with open(path + fn + ".txt", "w") as f:
        for (q, ints) in zip(q, ints):
            f.write("{0} {1}\n".format(q, ints))
    f.close()    
    
#-----------------------------------------------------------------------------
def importSim(path, fn):
    '''
    Imports text file with simulated XRD data

    Parameters
    ----------
    path : str
        File directory
    fn : str
        File name

    Returns
    -------
    q : array
        Calculated Q values
    ints : array
        Calculated intensity values

    '''
    q, ints = np.loadtxt(path + fn + ".txt", unpack=True, dtype=float)
    return q, ints

#-----------------------------------------------------------------------------
def loadCif(path, fn):
    '''
    Imports cif file as diffpy Structure object

    Parameters
    ----------
    path : str
        File directory
    fn : str
        File name

    Returns
    -------
    struct : Structure
        Structure object

    '''
    struct = df.Crystal(path + fn + ".cif")
    return struct

#-----------------------------------------------------------------------------
def norm(ints):
    '''
    Normalize intensity values

    Parameters
    ----------
    ints : list (float)
        Intensity values

    Returns
    -------
    norm_ints : array
        Normalized intensity values

    '''
    norm_ints = (ints / np.max(ints))
    return norm_ints


#-----------------------------------------------------------------------------
def hklSim(path, cif, wl, tt_max):
    '''
    Calculates (hkl) reflections

    Parameters
    ----------
    path : str
        File directory where cif is store
    cif : str
        File name
    wl : float
        Instrument wavelength (A)
    tt_max : float
        Maximum 2theta (deg)

    Returns
    -------
    reflections : list (float)
        Calculated (hkl) reflections, 2theta, and intensity values

    '''
    # load cif
    struct = df.Crystal(path + cif + ".cif")
    
    # setup diffraction parameters
    e_kev = df.fc.wave2energy(wl)
    struct.Scatter.setup_scatter(scattering_type="xray", energy_kev=e_kev, 
                                 min_twotheta=0, max_twotheta=tt_max)
    
    # generate (hkl) reflections
    reflections = struct.Scatter.print_all_reflections(min_intensity=0.1)
    
    return reflections


#------------------------------------------------------------------------------
def savehkl(path, fn, data):
    '''
    Saves (hkl) data as text file

    Parameters
    ----------
    path : str
        File directory
    fn : str
        File name
    data : list (float)
        (hkl) data

    Returns
    -------
    None

    '''
    with open(path + fn + "_hkl.txt", "w") as f:
        f.write(data)
    f.close()
    
#------------------------------------------------------------------------------
def diffCurve(q1, q2, ints1, ints2): 
    '''
    Calculates difference curve between two datasets

    Parameters
    ----------
    q1 : array
        Q data from dataset one
    q2 : array
        Q data from dataset two
    ints1 : array
        Intensity data from dataset one
    ints2 : array
        Intensity data from dataset two

    Returns
    -------
    diff_q : array
        Q data of difference curve
    diff_ints : array
        Intensity data of difference curve

    '''
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
def saveDiffCurve(q, ints, path, fn):
    '''
    Saves difference curve to text file

    Parameters
    ----------
    q : array
        Q data of difference curve
    ints : array
        Intensity data of difference curve
    path : str
        File directory
    fn : str
        File name

    Returns
    -------
    None

    '''
    with open(path + fn + ".txt", "w") as f:
        for (q, ints) in zip(q, ints):
            f.write("{0} {1}\n".format(q, ints))
    f.close()
    
    
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' XRD SIMULATION PLOTTING FUNCTIONS '''
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def tt_to_q(twotheta, wavelength):
    Q = 4 * np.pi * np.sin((twotheta * np.pi)/360) / wavelength
    return Q

#------------------------------------------------------------------------------
def saveFig(plot, path, fn):
    '''
    Save figure as a .png file

    Parameters
    ----------
    plot : figure
        Name of plot
    path : str
        File directory path
    fn : str
        File name to save to

    Returns
    -------
    None

    '''
    save = path + fn + ".png"
    plot.savefig(save, bbox_inches="tight", pad_inches=0.2, dpi=1000)

#------------------------------------------------------------------------------
def gradientGen(start_hex, end_hex, num):
    '''
    Generates color gradient

    Parameters
    ----------
    start_hex : str
        Hex code for first gradient color, format "#000000"
    end_hex : str
        Hex code for final gradient color, format "#000000"
    num : int
        Number of colors to generate

    Returns
    -------
    colors_list : list (Color)
        Hex codes, will need to use .hex() to retrieve as string

    '''
    start_color = Color(start_hex)
    end_color = Color(end_hex)
    
    colors_list = list(start_color.range_to(end_color, num))
    
    return colors_list

#------------------------------------------------------------------------------
def gradientGen2D(top_lf, top_rt, bott_lf, bott_rt, rows, cols):
    '''
    Generates a 2D color gradient

    Parameters
    ----------
    top_lf : str
        Hex code for top left gradient color, format "#000000"
    top_rt : str
        Hex code for top right gradient color, format "#000000"
    bott_lf : str
        Hex code for bottom left gradient color, format "#000000"
    bott_rt : str
        Hex code for bottom right gradient color, format "#000000"
    rows : int
        Number of rows
    cols : int
        Number of columns.

    Returns
    -------
    color_list : list (Color)
        Hex codes in 2D array format, will need to use .hex() to retrieve as string

    '''
    # get colors from hex codes
    tl = Color(top_lf)
    tr = Color(top_rt)
    bl = Color(bott_lf)
    br = Color(bott_rt)
    
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
def simStack(expt_q, expt_ints, num, q_list, ints_list, x_lim, y_lim, wl,
                    start_hex, end_hex, labels=None, label_offsets=None):
    '''
    Generates stacked plot of simulated vs. experimental data

    Parameters
    ----------
    expt_q : array
        Experimental Q data
    expt_ints : array
        Experimental intensity data (normalized)
    num : int
        Total number of datasets
    q_list : list (array)
        Arrays of simulated Q datasets
    ints_list : list (array)
        Arrays of simulated intensity datasets (normalized)
    x_lim : list (float)
        Tuple with x-axis minimum and maximum
    y_lim : list (float)
        Tuple with y-axis minimum and maximum
    wl : float
        Instrument wavelength (A)
    start_hex : str
        Hex code for initial gradient color, format "#000000"
    end_hex : str
        Hex code for final gradient color, format "#000000"
    labels : list (str), optional
        Labels for each dataset. The default is None.
    label_offsets : list (float), optional
        Tuple with offsets from x-axis maximum and vertical spacing from data
        for text labels. The default is None.

    Returns
    -------
    None

    '''
    
    fig, (p) = plt.subplots(num, 1, figsize=(8, 8))
    
    expt_min = np.min(expt_ints)
    
    # generate color gradient
    if start_hex == False:
        start_hex = "#00C6BF"
    if end_hex == False:
        end_hex = "#B430C2"
    g = gradientGen(start_hex, end_hex, num)
    
    # plot data
    for i in range(0, num):
        p[i].scatter(expt_q, expt_ints, color="black", label="Observed", marker=".", s=8)
    for i in range(0, num):
        p[i].plot(q_list[i], ints_list[i]+expt_min, color=g[i].hex, linewidth="2.5")
        
    # set axis limits
    for i in range(num):
        p[i].set_xlim(x_lim)
        p[i].set_ylim(y_lim)
        p[i].tick_params(axis="both", labelsize="14")
        if i < (num-1):
            p[i].get_xaxis().set_visible(False)
    
    # set axis labels
    x_label = r"Q (\AA" r"$^{-1}$, $\lambda=$" + str(wl) + r" \AA)"
    y_label = "Intensity (counts, normalized)"
    
    p[-1].set_xlabel(x_label, fontsize=16)
    fig.supylabel(y_label, fontsize=16)
    
    # add stack labels 
    if labels is not None:
        for i in range(num):
            p[i].text(x_lim[1] - label_offsets[0], label_offsets[1],
                    labels[i], color=g[i].hex, fontsize="16", ha="right", va="top")
    
    # add legend
    p[0].legend(handlelength=0.25, fontsize="14")
    
    plt.subplots_adjust(hspace=0.05) 
    
    return(p)

#------------------------------------------------------------------------------
def diffStack(num, q_list, ints_list, x_lim, y_lim, wl, start_hex, end_hex, 
               labels=None, label_offsets=None):
    '''
    Generates stacked plot of sim vs. expt difference curves

    Parameters
    ----------
    num : int
        Total number of datasets
    q_list : list (array)
        Arrays of difference curve Q datasets
    ints_list : list (array)
        Arrays of difference curve intensity datasets (normalized)
    x_lim : list (float)
        Tuple with x-axis minimum and maximum
    y_lim : list (float)
        Tuple with y-axis minimum and maximum
    wl : float
        Instrument wavelength (A)
    start_hex : str
        Hex code for initial gradient color, format "#000000"
    end_hex : str
        Hex code for final gradient color, format "#000000"
    labels : list (str), optional
        Labels for each dataset. The default is None.
    label_offsets : list (float), optional
        Tuple with offsets from x-axis maximum and vertical spacing from data
        for text labels. The default is None.

    Returns
    -------
    None

    '''
    
    fig, (p) = plt.subplots(num, 1, figsize=(8, 8))
    
    # generate color gradient
    if start_hex == False:
        start_hex = "#00C6BF"
    if end_hex == False:
        end_hex = "#B430C2"
    g = gradientGen(start_hex, end_hex, num)
    
    # plot data
    for i in range(0, num):
        p[i].plot(q_list[i], ints_list[i], color="#BEBEBE", linewidth="2")
        
    # set axis limits
    for i in range(num):
        p[i].set_xlim(x_lim)
        p[i].set_ylim(y_lim)
        p[i].tick_params(axis="both", labelsize="14")
        if i < (num-1):
            p[i].get_xaxis().set_visible(False)
    
    # set axis labels
    x_label = r"Q (\AA" r"$^{-1}$, $\lambda=$" + str(wl) + r" \AA)"
    y_label = "Intensity (counts, normalized)"
    
    p[-1].set_xlabel(x_label, fontsize=16)
    fig.supylabel(y_label, fontsize=16)
    
    # add stack labels 
    if labels is not None:
        for i in range(num):
            p[i].text(x_lim[1] - label_offsets[0], label_offsets[1],
                    labels[i], color=g[i].hex, fontsize="16", ha="right", va="top")
    
    plt.subplots_adjust(hspace=0.05) 
    
    return(p)

#------------------------------------------------------------------------------
def compare_uf_flt(expt_q, expt_ints, uf_q, uf_ints, flt_q, flt_ints, x_lim, 
                   y_lim, wl, calcDiff=False, diff_q=None, diff_ints=None):
    '''
    Generates plot of unfaulted sim vs. expt and one faulted sim vs. expt

    Parameters
    ----------
    expt_q : array
        Experimental Q data
    expt_ints : array
        Experimental intensity data (normalized)
    uf_q : array
        Unfaulted supercell Q data
    uf_ints : array
        Unfaulted supercell intensity data
    flt_q : array
        Faulted supercell Q data
    flt_ints : array
        Faulted supercell intensity data
    x_lim : list (float)
        Tuple with x-axis minimum and maximum
    y_lim : list (float)
        Tuple with y-axis minimum and maximum
    wl : float
        Instrument wavelength (A)

    Returns
    -------
    None

    '''
    
    fig, (p) = plt.subplots(1, 2, sharey=True, figsize=(14,7))
    
    c = ["#FA008E", "#2FF8B9"]
    expt_min = np.min(expt_ints)
    
    # calculate difference curves
    if calcDiff == True:
        uf_diff_q, uf_diff_ints = diffCurve(expt_q, uf_q, expt_ints, uf_ints + expt_min)
        flt_diff_q, flt_diff_ints = diffCurve(expt_q, flt_q, expt_ints, flt_ints + expt_min)
        
        diff_q = [uf_diff_q, flt_diff_q]
        diff_ints = [uf_diff_ints, flt_diff_ints]
    
    # plot expt data
    for i in range(2):
        p[i].scatter(expt_q, expt_ints, color="black", label="Observed", marker=".")
    
    # plot unfaulted data
    p[0].plot(uf_q, uf_ints, color=c[0], label="Unfaulted", linewidth="2")
    p[0].plot(diff_q[0], diff_ints[0]-np.max(diff_ints[0]), color="#BEBEBE", 
              label="Difference", linewidth="1")

    # plot faulted data
    p[1].plot(flt_q, flt_ints, color=c[1], label="Faulted", linewidth="2")
    p[1].plot(diff_q[1], diff_ints[1]-np.max(diff_ints[1]), color="#BEBEBE", 
              label="Difference", linewidth="1")
    
    # set axis limits
    for i in range(2):
        p[i].set_xlim(x_lim)
        p[i].set_ylim(y_lim)
        p[i].tick_params(axis="both", labelsize="14")
    
    p[1].get_yaxis().set_visible(False)

    # set axis labels
    x_label = r"Q (\AA" r"$^{-1}$, $\lambda=$" + str(wl) + r" \AA)"
    y_label = "Intensity (counts, normalized)"
    
    fig.supxlabel(x_label, fontsize=16)
    p[0].set_ylabel(y_label, fontsize=16)
    
    for i in range(2):
        p[i].legend(handlelength=1, fontsize="14", ncols=3, loc="upper center")
    
    plt.subplots_adjust(wspace=0.05)
    
    return(p)
    
#------------------------------------------------------------------------------
def addFltParamBox(plot, n_stacks, pos, ha="left", va="top", p=None, s=None, 
                      color=None):
    '''
    Adds a text box with supercell parameters

    Parameters
    ----------
    plot : Figure
        Plot to add text to
    n_stacks : int
        Number of stacks in supercell
    pos : list (float)
        Tuple with (x, y) position
    ha : str, optional
        Horizontal text alignment. The default is "left".
    va : str, optional
        Vertical text alignment. The default is "top".
    p : float, optional
        Stacking probability. The default is None.
    s : list (float), optional
        Stacking vector [x, y]. The default is None.
    color : str, optional
        Hex code for text box fill color, format "#000000". The default is None.

    Returns
    -------
    None

    '''
    
    if p is None:
        title = "Unfaulted"
    if p is not None:
        title = "Faulted"
        
    n_txt = re.sub("x", str(n_stacks), r"$N = x$")
    
    if p is not None:
        p_txt = re.sub("x", str(int(p*100)), r"$P = x \%$")
        
    if s is not None:
        s[0] = str(s[0])
        s[1] = str(s[1])
        
        isFrac = [False, False]
        
        for i in range(len(s)):
            if "/" in s[i]:
                isFrac[i] = True
                
        if isFrac[0] == False and isFrac[1] == False:
            vec_txt = r"$\vec{S} = \left[ x, y \right]$"
            var_list = ["x", "y"]
            
            for i in range(len(s)):
                sub_txt = re.sub(var_list[i], str(s[i]), vec_txt)
                vec_txt = sub_txt
            
        if isFrac[0] == False and isFrac[1] == True:
            vec_txt = r"$\vec{S} = \left[ x, \frac{y1}{y2} \right]$"
            var_list = ["x", "y1", "y2"]
            
            split_s = s[1].split("/")
            new_s = [s[0], split_s[0], split_s[1]]
            
            for i in range(len(new_s)):
                sub_txt = re.sub(var_list[i], str(new_s[i]), vec_txt)
                vec_txt = sub_txt
            
        if isFrac[0] == True and isFrac[1] == False:
            vec_txt = r"$\vec{S} = \left[ \frac{x1}{x2}, y \right]$"
            var_list = ["x1", "x2", "y"]
            
            split_s = s[0].split("/")
            new_s = [split_s[0], split_s[1], s[1]]
            
            for i in range(len(new_s)):
                sub_txt = re.sub(var_list[i], str(new_s[i]), vec_txt)
                vec_txt = sub_txt
            
        if isFrac[0] == True and isFrac[1] == True:
            vec_txt = r"$\vec{S} = \left[ \frac{x1}{x2}, \frac{y1}{y2} \right]$"
            var_list = ["x1", "x2", "y1", "y2"]
            
            split_sx = s[0].split("/")
            split_sy = s[1].split("/")
            new_s = [split_sx[0], split_sx[1], split_sy[0], split_sy[1]]
            
            for i in range(len(new_s)):
                sub_txt = re.sub(var_list[i], str(new_s[i]), vec_txt)
                vec_txt = sub_txt
        
    if p is None:
        txt = "\n".join((title, n_txt))
    if p is not None:
        txt = "\n".join((title, n_txt, p_txt, vec_txt))
    
    if color is None:
        color = "white"
    props = dict(boxstyle="round", facecolor=color, alpha=0.5, pad=0.3)
    
    plot.text(pos[0], pos[1], txt, bbox=props, ha=ha, va=va, fontsize="16", linespacing=2)
    
#------------------------------------------------------------------------------
def addPeakLabels(plot, hkl, x_pos, y_pos, color=None, size="14"):
    '''
    Adds (hkl) text labels to reflections

    Parameters
    ----------
    plot : Figure
        Plot to add text to
    hkl : list (str)
        (hkl) text labels
    x_pos : list (float)
        x-axis positions for labels
    y_pos : list (float)
        y-axis positions for labels
    color : str, optional
        Hex code for text color, format "#000000". The default is None.
    size : str, optional
        Font size. The default is "14".

    Returns
    -------
    None

    '''
    if color is None:
        color = "black"
    
    for i in range(len(hkl)):
        plot.text(x_pos[i], y_pos[i], hkl[i], color=color, ha="center", 
                  va="center", fontsize=size)

#------------------------------------------------------------------------------
def fitCompare(rows, cols, diff_q, diff_ints, x_lims, y_lim, wl, row_labels, 
                col_labels, row_label_adj, col_label_adj):
    '''
    Compare goodness of fit between datasets with a difference of difference curve

    Parameters
    ----------
    rows : int
        Number of rows
    cols : int
        Number of columns
    diff_q : list (float)
        List of Q datasets of difference of difference curves
    diff_ints : list (float)
        List of intensity datasets of difference of difference curves
    x_lims : list (float)
        List of tuples with x-axis minimums and maximums for each column
    y_lim : list (float)
        Tuple with y-axis minimum and maximum
    wl : float
        Instrument wavelength (A)
    row_labels : list (str)
        Text labels for rows
    col_labels : list (str)
        Text labels for columns
    row_label_adj : float
        Value to shift row labels along x-axis
    col_label_adj : float
        Value to shift column labels along y-axis

    Returns
    -------
    None

    '''
    
    fig, (p) = plt.subplots(rows, cols, figsize=(rows*2, cols))
    
    g = gradientGen2D("#00C6BF", "#009AE1", "#5D7AD3", "#B430C2", rows, cols)
    
    for row in range(rows):
        for col in range(cols):
            # plot data
            p[row][col].plot(diff_q[row], diff_ints[row], color=g[row][col].hex)
            
    # set axis limits
    for row in range(rows):
        for col in range(cols):
            p[row][col].set_xlim(x_lims[col])
            p[row][col].set_ylim(y_lim)
            p[row][col].tick_params(axis="both", labelsize="14")
            
    # format axes
    for row in range(rows - 1):
        for col in range(cols):
            p[row][col].get_xaxis().set_visible(False)
    
    for row in range(rows):
        for col in range(1, cols):
            p[row][col].get_yaxis().set_visible(False)
      
    for col in range(cols):
        p[-1][col].xaxis.set_major_locator(ticker.MultipleLocator(0.01))
    
    # set axis labels
    x_label = r"Q (\AA" r"$^{-1}$, $\lambda=$" + str(wl) + r" \AA)"
    y_label = r"Diff$_{\mathrm{UF}} -$ Diff$_{\mathrm{F}}$ (counts, normalized)"
    
    fig.supxlabel(x_label, fontsize=16)
    fig.supylabel(y_label, fontsize=16)
    
    # set plot labels
    y_mid = ((y_lim[1] - y_lim[0]) / 2) + y_lim[0]
    x_end = x_lims[-1][1]
    
    for row in range(rows):
        p[row][-1].text(x_end + row_label_adj, y_mid, row_labels[row],
                        color=g[row][-1].hex, fontsize="16", ha="left", va="center")
        
    for col in range(cols):
        x_mid = ((x_lims[col][1] - x_lims[col][0]) / 2) + x_lims[col][0]
        p[0][col].text(x_mid, y_lim[1] + col_label_adj, col_labels[col], color=g[0][col].hex, 
                       fontsize="16", ha="center", va="bottom")
        
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    
    return(p)





