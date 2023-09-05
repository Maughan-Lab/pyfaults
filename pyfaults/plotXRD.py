#########################################################################################
# pyfaults.plotXRD
# Author: Sinclair R. Combs
#########################################################################################

from colour import Color
import numpy as np
import re

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
from matplotlib.pyplot import rc
rc("text", usetex=True)
rc("font", **{"family":"sans-serif","sans-serif":["Helvetica"]},size="14")
rc("text.latex",preamble=r"\usepackage{sfmath}")

#########################################################################################
''' XRD plotting functions '''
#########################################################################################

''' Saves figure as PNG file '''
#----------------------------------------------------------------------------------------
def saveFig(plot, path, fn):
    '''
    Parameters
    ----------
    plot : figure
        Name of plot
    path : str
        File directory path
    fn : str
        File name to save to
    '''
    save = path + fn + ".png"
    plot.savefig(save, bbox_inches="tight", pad_inches=0.2, dpi=1000)
 
    
''' Normalizes intensity values '''
#----------------------------------------------------------------------------------------
def norm(ints):
    '''
    Parameters
    ----------
    ints : nparray
        Intensity values

    Returns
    -------
    norm_ints : nparray
        Normalized intensity values
    '''
    normInts = (ints / np.max(ints))
    return normInts
    

''' Generates color gradient '''
#----------------------------------------------------------------------------------------
def gradientGen(startHex, endHex, num):
    '''
    Parameters
    ----------
    startHex : str
        Hex code for first gradient color, format "#000000"
    endHex : str
        Hex code for final gradient color, format "#000000"
    num : int
        Number of colors to generate

    Returns
    -------
    colors_list : list (Color)
        Hex codes, will need to use .hex() to retrieve as string
    '''
    startColor = Color(startHex)
    endColor = Color(endHex)
    
    colorsList = list(startColor.range_to(endColor, num))
    
    return colorsList


''' Generates 2D color gradient '''
#----------------------------------------------------------------------------------------
def gradientGen2D(cornerColors, rows, cols):
    '''
    Parameters
    ----------
    cornerColors : list (str)
        Hex codes for colors in each corner of gradient map, format "#000000".
        List colors in position order: top left, top right, bottom left, bottom
        right.
    rows : int
        Number of rows
    cols : int
        Number of columns.

    Returns
    -------
    color_list : list (Color)
        Hex codes in 2D array format, will need to use .hex() to retrieve as string
    '''
    tl = Color(cornerColors[0])
    tr = Color(cornerColors[1])
    bl = Color(cornerColors[2])
    br = Color(cornerColors[3])
    
    # generate gradient from upper left to lower left
    startColColors = list(tl.range_to(bl, rows))
    
    # generate gradient from upper right to lower lright
    endColColors = list(tr.range_to(br, rows))
    
    # generate remaining gradient colors
    colorList = []
    for i in range(rows):
        c1 = startColColors[i]
        c2 = endColColors[i]
        currRowColors = list(c1.range_to(c2, cols))
        colorList.append(currRowColors)
    
    return colorList


''' Generates stacked plot with experimental data '''
#----------------------------------------------------------------------------------------
def exptStackedPlot(num, expt, qVals, intsVals, x_lim, y_lim, wl, gradient=None, 
             labels=None, labelOffsets=None, normalized=False):
    '''
    Parameters
    ----------
    num : int
        Total number of datasets
    expt : list (nparray)
        Tuple of arrays of experimental Q data and intensity data
    qVals : list (nparray)
        List of Q datasets
    intsVals : list (nparray)
        List of intensity datasets
    x_lim : list (float)
        Tuple with x-axis minimum and maximum
    y_lim : list (float)
        Tuple with y-axis minimum and maximum
    wl : float
        Instrument wavelength (A)
    gradient : list (str), optional
        Tuple of start and end values of color gradient, format "#000000". The 
        default is ("#00C6BF", "#B430C2").
    labels : list (str), optional
        Labels for each dataset. The default is None.
    labelOffsets : list (float), optional
        Tuple with offsets from x-axis maximum and vertical spacing from data
        for text labels. The default is None.
    normalized : bool, optional
        Set to True if intensity is normalized. The default is False.
    '''
    
    fig, (p) = plt.subplots(num, 1, figsize=(8, 8))
    
    exptMin = np.min(expt[1])
    
    # generate gradient
    if gradient is not None:
        g = gradientGen(gradient[0], gradient[1], num)
    elif gradient is None:
        g = gradientGen("#00C6BF", "#B430C2", num)
    
    # plot data
    for i in range(num):
        p[i].scatter(expt[0], expt[1], color="black", label="Observed", marker=".", s=8)
    for i in range(num):
        p[i].plot(qVals[i], intsVals[i]+exptMin, color=g[i].hex, linewidth="2.5")
        
    # set axis limits
    for i in range(num):
        p[i].set_xlim(x_lim)
        p[i].set_ylim(y_lim)
        p[i].tick_params(axis="both", labelsize="14")
        if i < (num-1):
            p[i].get_xaxis().set_visible(False)
    
    # set axis labels
    x_label = r"Q (\AA" r"$^{-1}$, $\lambda=$" + str(wl) + r" \AA)"
    if normalized == True:
        y_label = "Intensity (counts, normalized)"
    elif normalized == False:
        y_label = "Intensity (counts)"
    
    p[-1].set_xlabel(x_label, fontsize=16)
    fig.supylabel(y_label, fontsize=16)
    
    # add stack labels 
    if labels is not None:
        for i in range(num):
            p[i].text(x_lim[1] - labelOffsets[0], labelOffsets[1],
                    labels[i], color=g[i].hex, fontsize="16", ha="right", va="top")
    
    # add legend
    p[0].legend(handlelength=0.25, fontsize="14")
    
    plt.subplots_adjust(hspace=0.05) 
    
    return(p)


''' Generates stacked plot '''
#----------------------------------------------------------------------------------------
def stackedPlot(num, qVals, intsVals, x_lim, y_lim, wl, gradient=None, 
                labels=None, labelOffsets=None, normalized=False):
    '''
    Parameters
    ----------
    num : int
        Total number of datasets
    qVals : list (nparray)
        List of Q datasets
    intsVals : list (nparray)
        List of intensity datasets
    x_lim : list (float)
        Tuple with x-axis minimum and maximum
    y_lim : list (float)
        Tuple with y-axis minimum and maximum
    wl : float
        Instrument wavelength (A)
    gradient : list (str), optional
        Tuple of start and end values of color gradient, format "#000000". The 
        default is ("#00C6BF", "#B430C2").
    labels : list (str), optional
        Labels for each dataset. The default is None.
    labelOffsets : list (float), optional
        Tuple with offsets from x-axis maximum and vertical spacing from data
        for text labels. The default is None.
    normalized : bool, optional
        Set to True if intensity is normalized. The default is False.
    '''
    
    fig, (p) = plt.subplots(num, 1, figsize=(8, 8))
    
    # generate gradient
    if gradient is not None:
        g = gradientGen(gradient[0], gradient[1], num)
    elif gradient is None:
        g = gradientGen("#00C6BF", "#B430C2", num)

    # plot data
    for i in range(num):
        p[i].plot(qVals[i], intsVals[i], color=g[i].hex, linewidth="2.5")
        
    # set axis limits
    for i in range(num):
        p[i].set_xlim(x_lim)
        p[i].set_ylim(y_lim)
        p[i].tick_params(axis="both", labelsize="14")
        if i < (num-1):
            p[i].get_xaxis().set_visible(False)
    
    # set axis labels
    x_label = r"Q (\AA" r"$^{-1}$, $\lambda=$" + str(wl) + r" \AA)"
    if normalized == True:
        y_label = "Intensity (counts, normalized)"
    elif normalized == False:
        y_label = "Intensity (counts)"
    
    p[-1].set_xlabel(x_label, fontsize=16)
    fig.supylabel(y_label, fontsize=16)
    
    # add stack labels 
    if labels is not None:
        for i in range(num):
            p[i].text(x_lim[1] - labelOffsets[0], labelOffsets[1],
                    labels[i], color=g[i].hex, fontsize="16", ha="right", va="top")
    
    plt.subplots_adjust(hspace=0.05) 
    
    return(p)


''' Generates dual plot. Plot 1 shows experimental XRD and unfaulted supercell
XRD. Plot 2 shows experimental XRD and faulted supercell XRD. '''
#----------------------------------------------------------------------------------------
def compareUFtoFLT(expt, UF, FLT, UFdiff, FLTdiff, nStacks, x_lim, y_lim, wl, 
                   prob, sVec, boxAdj=None, legendAdj=None, colors=None, 
                   normalized=False, diffOffset=0):
    '''
    Parameters
    ----------
    expt : list (nparray)
        Tuple of arrays of experimental Q data and intensity data
    UF : list (nparray)
        Tuple of arrays of unfaulted supercell Q data and intensity data
    FLT : list (nparray)
        Tuple of arrays of faulted supercell Q data and intensity data
    UFdiff : list (nparray)
        Tuple of arrays of expt vs UF difference Q data and intensity data
    FLTdiff : list (nparray)
        Tuple of arrays of expt vs FLT difference Q data and intensity data
    nStacks : int
        Number of stacks in supercells
    x_lim : list (float)
        Tuple with x-axis minimum and maximum
    y_lim : list (float)
        Tuple with y-axis minimum and maximum
    wl : float
        Instrument wavelength (A)
    prob : float
        Stacking probability
    sVec : list (str)
        List of stacking vector components as strings
    boxAdj : list (float), optional
        Tuple of parameter box position displacement from upper left corner.
    legendAdj : list (float), optional
        Tuple of legend position displacement. Default is (0.01, 0.01).
    colors : list (str), optional
        Tuple of color hex codes for unfaulted and faulted curves, format 
        "#000000". The default is ("#00C6BF", "#B430C2").
    normalized : bool, optional
        Set to True if intensity is normalized. The default is False.
    diffOffset : float, optional
        Additional adjustment parameter for difference curve placement. The
        default is 0.
    '''
    
    fig, (p) = plt.subplots(1, 2, sharey=True, figsize=(14,7))
    
    if colors is not None:
        c = colors
    if colors is None:
        c = ("#00C6BF", "#B430C2")
        
    exptMin = np.min(expt[1])
    
    # plot expt data
    for i in range(2):
        p[i].scatter(expt[0], expt[1], color="black", label="Observed", marker=".")
    
    # plot unfaulted data
    p[0].plot(UF[0], UF[1] + exptMin, color=c[0], label="Unfaulted", linewidth="3")
    p[0].plot(UFdiff[0], UFdiff[1] + diffOffset, color="#BEBEBE", 
              label="Difference", linewidth="1")

    # plot faulted data
    p[1].plot(FLT[0], FLT[1] + exptMin, color=c[1], label="Faulted", linewidth="3")
    p[1].plot(FLTdiff[0], FLTdiff[1] + diffOffset, color="#BEBEBE", 
              label="Difference", linewidth="1")
    
    # set axis limits
    for i in range(2):
        p[i].set_xlim(x_lim)
        p[i].set_ylim(y_lim)
        p[i].tick_params(axis="both", labelsize="16")
    
    p[1].get_yaxis().set_visible(False)

    # set axis labels
    x_label = r"Q (\AA" r"$^{-1}$, $\lambda=$" + str(wl) + r" \AA)"
    if normalized == True:
        y_label = "Intensity (counts, normalized)"
    elif normalized == False:
        y_label = "Intensity (counts)"
    
    fig.supxlabel(x_label, fontsize=18)
    p[0].set_ylabel(y_label, fontsize=18)
    
    # set legend handles
    obsHandle = mlines.Line2D([], [], color="white", label="Observed", marker=".", 
                     mfc="black", ms=15)
    UFlabel = "\n".join(("Unfaulted", re.sub("x", str(nStacks), r"$N = x$")))
    UFhandle = mlines.Line2D([], [], color=c[0], label=UFlabel)
    FLTlabel = "\n".join(("Faulted", re.sub("x", str(nStacks), r"$N = x$")))
    FLThandle = mlines.Line2D([], [], color=c[1], label=FLTlabel)
    diffHandle = mlines.Line2D([], [], color="#BEBEBE", label="Difference")
    
    if legendAdj is not None:
        la = legendAdj
    elif legendAdj is None:
        la = (0.01, 0.01)
    
    # add legend
    p[1].legend(handles=[obsHandle, UFhandle, FLThandle, diffHandle], 
                handlelength=1, fontsize="16", loc="upper left", 
                bbox_to_anchor=(1 + la[0], 1 + la[1]))
    
    # add fault parameters box
    probText = re.sub("x", str(int(prob*100)), r"$P = x \%$")
    sVecText = r"$\vec{S} = \left[ x, y, z \right]$"
    varList = ["x", "y", "z"]
    for i in range(3):
        subText = re.sub(varList[i], sVec[i], sVecText)
        sVecText = subText
    
    if boxAdj is not None:
        ba = boxAdj
    elif boxAdj is None:
        ba = (0.01, -0.01)
        
    params = "\n".join((sVecText, probText))
    props = dict(boxstyle="round", facecolor="white", alpha=0.5, pad=0.3)
    p[1].text(x_lim[0] + ba[0], y_lim[1] + ba[1], params, bbox=props, 
              ha="left", va="top", fontsize="16", linespacing=2)
    
    plt.subplots_adjust(wspace=0.05)
    
    return(p)

''' Add labels to specific (hkl) reflections '''
#----------------------------------------------------------------------------------------
def addPeakLabels(plot, hkl, x_pos, y_pos, color=None, size="14"):
    '''
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
    '''
    if color is None:
        color = "black"
    
    for i in range(len(hkl)):
        plot.text(x_pos[i], y_pos[i], hkl[i], color=color, ha="center", 
                  va="center", fontsize=size)


''' Compare goodness of fit between data sets using difference of differences.
Generates plot with rows corresponding to different models and columns
corresponding to different reflections of interest.'''
#----------------------------------------------------------------------------------------
def fitCompare(rows, cols, diffQ, diffInts, x_lims, y_lim, wl, rowLabels, 
                colLabels, rowLabelAdj=None, colLabelAdj=None, xLabelAdj=None, 
                yLabelAdj=None, gradient=None, normalized=False):
    '''
    Parameters
    ----------
    rows : int
        Number of rows
    cols : int
        Number of columns
    diffQ : list (nparray)
        Lists of Q datasets of difference of difference curves, format as list
        of lists where each row entry has [c1, c2, c3, ...]
    diffInts : list (nparray)
        List of intensity datasets of difference of difference curves, format
        as list of lists where each row entry has [c1, c2, c3, ...]
    x_lims : list (float)
        List of tuples with x-axis minimums and maximums for each column
    y_lim : list (float)
        Tuple with y-axis minimum and maximum
    wl : float
        Instrument wavelength (A)
    rowLabels : list (str)
        Text labels for rows
    colLabels : list (str)
        Text labels for columns
    rowLabelAdj : float, optional
        Adjustment paramter for row label position along x-axis. The default 
        is 0.01.
    colLabelAdj : float, optional
        Adjustment paramter for column label position along y-axis. The default
        is 0.01.
    xLabelAdj : float, optional
        Adjustment parameter for x-axis label position. The default is 0.
    yLabelAdj : float, optional
        Adjustment parameter for y-axis label position. The default is 0.
    gradient : list (str)
        Hex codes for colors in each corner of gradient map, format "#000000".
        List colors in position order: top left, top right, bottom left, bottom
        right. The default is ["#00C6BF", "#009AE1", "#5D7AD3", "#B430C2"].
    normalized : bool, optional
        Set to True if intensity is normalized. The default is False.
    '''
    
    fig, (p) = plt.subplots(rows, cols, figsize=(rows*2, cols))
    
    if gradient is not None:
        g = gradientGen2D(gradient, rows, cols)
    elif gradient is None:
        g = gradientGen2D(["#00C6BF", "#009AE1", "#5D7AD3", "#B430C2"], rows, cols)
        
    if rows == 1:
        for col in range(cols):
            p[col].plot(diffQ[0][col], diffInts[0][col], color=g[0][col].hex)
    elif cols == 1:
        for row in range(rows):
            p[row].plot(diffQ[row][0], diffInts[row][0], color=g[row][0].hex)
    else:
        for row in range(rows):
            for col in range(cols):
                # plot data
                p[row][col].plot(diffQ[row][col], diffInts[row][col], 
                                 color=g[row][col].hex)
            
    # set axis limits
    for row in range(rows):
        for col in range(cols):
            p[row][col].set_xlim(x_lims[col][0], x_lims[col][1])
            p[row][col].set_ylim(y_lim)
            p[row][col].tick_params(axis="both", labelsize="14")
            
    # format axes
    for row in range(rows - 1):
        for col in range(cols):
            p[row][col].get_xaxis().set_visible(False)
    
    for row in range(rows):
        for col in range(1, cols):
            p[row][col].get_yaxis().set_visible(False)
    
    # set axis labels
    x_label = r"Q (\AA" r"$^{-1}$, $\lambda=$" + str(wl) + r" \AA)"
    if normalized == True:
        y_label = r"Diff$_{\mathrm{UF}} -$ Diff$_{\mathrm{F}}$ (counts, normalized)"
    elif normalized == False:
        y_label = r"Diff$_{\mathrm{UF}} -$ Diff$_{\mathrm{F}}$ (counts)"
    
    fig.supxlabel(x_label, fontsize=16)
    fig.supylabel(y_label, fontsize=16)
    
    # set plot labels
    y_mid = ((y_lim[1] - y_lim[0]) / 2) + y_lim[0]
    x_end = x_lims[-1][1]
    
    if rowLabelAdj is None:
        rowLabelAdj = 0.01
    if colLabelAdj is None:
        colLabelAdj = 0.01
    
    for row in range(rows):
        p[row][-1].text(x_end + rowLabelAdj, y_mid, rowLabels[row],
                        color=g[row][-1].hex, fontsize="16", ha="left", va="center")
        
    for col in range(cols):
        x_mid = ((x_lims[col][1] - x_lims[col][0]) / 2) + x_lims[col][0]
        p[0][col].text(x_mid, y_lim[1] + colLabelAdj, colLabels[col], color=g[0][col].hex, 
                       fontsize="16", ha="center", va="bottom")
        
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    
    return(p)


''' Set spacing of tick marks on plot axis '''
#----------------------------------------------------------------------------------------
def setTickSpacing(plot, axis, interval):
    '''
    Parameters
    ----------
    plot : Figure
        Plot to adjust tick spacing on
    axis : str
        Axis to adjust spacing on; set to "x", "y", or "both"
    interval : float
        Spacing interval
    '''
    
    if axis == "x":
        plot.xaxis.set_major_locator(ticker.MultipleLocator(interval))
    elif axis == "y":
        plot.yaxis.set_major_locator(ticker.MultipleLocator(interval))
    elif axis == "both":
        plot.xaxis.set_major_locator(ticker.MultipleLocator(interval))
        plot.yaxis.set_major_locator(ticker.MultipleLocator(interval))
        
''' Adjust axis label position '''
#----------------------------------------------------------------------------------------


