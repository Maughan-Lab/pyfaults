##################################################################################
# Author: Sinclair R. Combs
##################################################################################

from colour import Color
import numpy as np
import re

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
from matplotlib.pyplot import rc
rc('text', usetex=True)
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']},size='14')
rc('text.latex',preamble=r'\usepackage{sfmath}')


''' XRD PLOTTING METHODS '''
#---------------------------------------------------------------------------------
# saves to PNG -------------------------------------------------------------------
#---------------------------------------------------------------------------------
def saveFig(plot, path, fn):
    '''
    Parameters
    ----------
    plot
        Figure : unique identifier for plot
    path
        str : file save directory
    fn
        str : file name
    '''
    save = path + fn + '.png'
    plot.savefig(save, bbox_inches='tight', pad_inches=0.2, dpi=1000)

#---------------------------------------------------------------------------------
# generates 1D color gradient ----------------------------------------------------
#---------------------------------------------------------------------------------
def gradientGen(startHex, endHex, num):
    '''
    Parameters
    ----------
    startHex
        str : hex code for first color, format '#000000'
    endHex
        str : hex code for last color, format '#000000'
    num
        int : number of colors in gradient

    Returns
    -------
    colorList
        list (Color) : gradient hex codes
        Note: must use .hex() to retrieve codes as strings
    '''
    startColor = Color(startHex)
    endColor = Color(endHex)
    colorList = list(startColor.range_to(endColor, num))
    return colorList

#---------------------------------------------------------------------------------
# generates 2D color gradient ----------------------------------------------------
#---------------------------------------------------------------------------------
def gradientGen2D(cornerColors, rows, cols):
    '''
    Parameters
    ----------
    cornerColors
        list (str) : hex codes for corners of gradient map, format '#000000'
    rows
        int : number of rows
    cols
        int : number of columns

    Returns
    -------
    colorList
        list (Color) : gradient hex codes
        Note: must use .hex() to retrieve codes as strings
    '''
    tl = Color(cornerColors[0])
    tr = Color(cornerColors[1])
    bl = Color(cornerColors[2])
    br = Color(cornerColors[3])
    
    startColColors = list(tl.range_to(bl, rows))
    endColColors = list(tr.range_to(br, rows))
    colorList = []
    for i in range(rows):
        c1 = startColColors[i]
        c2 = endColColors[i]
        currRowColors = list(c1.range_to(c2, cols))
        colorList.append(currRowColors)
    return colorList

#---------------------------------------------------------------------------------
# generates stacked plot ---------------------------------------------------------
#---------------------------------------------------------------------------------
def stackedPlot(num, qVals, intsVals, xLim, yLim, wl, expt=None,
                    gradient=None, labels=None, labelOffsets=(0,0), 
                    normalized=False):
    '''
    Parameters
    ----------
    num
        int : number of datasets
    qVals
        list (nparray) : list of Q value datasets
    intsVals
        list (nparray) : list of intensity value datasets
    xLim
        list (float) : x-axis minimum and maximum, tuple
    yLim
        list (float) : y-axis minimum and maximum, tuple
    wl
        float : instrument wavelength (A)
    expt
        list (nparray) [optional] : experimental Q and intensity values
    gradient
        list (str) [optional] : first and last gradient color hex codes
        The default is ['#00C6BF', '#B430C2']
    labels
        list (str) [optional] : list of text labels for each dataset
        The default is None
    labelOffsets
        list (float) [optional] : adjustment parameters for text labels
        The default is (0,0)
        Note: offsets from x-axis maximum and stacked datasets
    normalized
        bool [optional] : set to True if intensity values are normalized
        The default is False
    '''
    
    fig, (p) = plt.subplots(num, 1, figsize=(8, 8))
    exptMin = np.min(expt[1])
    
    # generate gradient
    if gradient is not None:
        g = gradientGen(gradient[0], gradient[1], num)
    elif gradient is None:
        g = gradientGen('#00C6BF', '#B430C2', num)
    
    # plot experimental data
    if expt is not None:
        for i in range(num):
            p[i].scatter(expt[0], expt[1], 
                         color='black', 
                         label='Observed', 
                         marker='.', s=8)
        
    # plot stacked data
    for i in range(num):
        p[i].plot(qVals[i], intsVals[i]+exptMin, 
                  color=g[i].hex, 
                  linewidth='2.5')
        # set axis limits
        p[i].set_xlim(xLim)
        p[i].set_ylim(yLim)
        p[i].tick_params(axis='both', labelsize='14')
        if i < (num-1):
            p[i].get_xaxis().set_visible(False)
    
    # set x-axis label
    x_label = r'Q (\AA" r"$^{-1}$, $\lambda=$' + str(wl) + r' \AA)'
    p[-1].set_xlabel(x_label, fontsize=16)
    
    # set y-axis label
    if normalized == True:
        y_label = 'Intensity (counts, normalized)'
    elif normalized == False:
        y_label = 'Intensity (counts)'
    fig.supylabel(y_label, fontsize=16)
    
    # add additional labels 
    if labels is not None:
        for i in range(num):
            p[i].text(xLim[1] - labelOffsets[0], labelOffsets[1],
                      labels[i], 
                      color=g[i].hex, 
                      fontsize='16', 
                      ha='right', va='top')
    
    # add legend
    p[0].legend(handlelength=0.25, fontsize='14')
    
    plt.subplots_adjust(hspace=0.05) 
    return(p)

#---------------------------------------------------------------------------------
# generates dual plot of a) experimental vs unfaulted supercell simulation XRD ---
# and b) experimental vs faulted supercell simulation XRD ------------------------
#---------------------------------------------------------------------------------
def compareUFtoFLT(expt, UF, FLT, UFdiff, FLTdiff, nStacks, prob, sVec, 
                   xLim, yLim, wl, simColors=None, normalized=False,
                   boxAdj=None, legendAdj=None, diffAdj=None):
    '''
    Parameters
    ----------
    expt
        list (nparray) : experimental Q and intensity values
    UF
        list (nparray) : unfaulted supercell Q and intensity values
    FLT
        list (nparray) : faulted supercell Q and intensity values
    UFdiff
        list (nparray) : experimental vs unfaulted difference curve values
    FLTdiff
        list (nparray) : experimental vs faulted difference curve values
    nStacks
        int : number of stacks in supercell
    prob
        float : stacking fault probability
    sVec
        list (str) : stacking vector [x, y, z] components in string format 
    xLim
        list (float) : x-axis minimum and maximum, tuple
    yLim
        list (float) : y-axis minimum and maximum, tuple
    wl
        float : instrument wavelength (A)
    simColors
        list (str) [optional] : color hex codes for simulated XRD patterns
        The default is ["#00C6BF", "#B430C2"]
    normalized
        bool [optional] : set to True if intensity values are normalized
        The default is False
    boxAdj
        list (float) [optional] : adjustment parameters for text box location
        The default is (0.01, -0.01)
        Note: offsets from upper left corner
    legendAdj
        list (float) [optional] : adjustment parameters for legend location
        The default is (0.01, 0)
    diffAdj
        float [optional] : adjustment parameter for difference curve location
        The default is 0
    '''
    
    fig, (p) = plt.subplots(1, 2, sharey=True, figsize=(14,7))
    exptMin = np.min(expt[1])
    
    if simColors is None:
        simColors = ('#00C6BF', '#B430C2')
    if boxAdj is None:
        boxAdj = (0.01, -0.01)
    if legendAdj is None:
        legendAdj = (0.01, 0)
    if diffAdj is None:
        diffAdj = 0
        
    
    # plot experimental data
    for i in range(2):
        p[i].scatter(expt[0], expt[1], 
                     color='black', 
                     label='Observed', 
                     marker='.')
        # set axis limits
        p[i].set_xlim(xLim)
        p[i].set_ylim(yLim)
        p[i].tick_params(axis='both', labelsize='16')
        p[1].get_yaxis().set_visible(False)
    
    # plot unfaulted supercell data
    p[0].plot(UF[0], UF[1] + exptMin, 
              color=simColors[0], 
              label='Unfaulted', 
              linewidth='3')
    # plot experimental vs unfaulted difference
    p[0].plot(UFdiff[0], UFdiff[1] + diffAdj, 
              color='#BEBEBE', 
              label='Difference', 
              linewidth='1')

    # plot faulted supercell data
    p[1].plot(FLT[0], FLT[1] + exptMin, 
              color=simColors[1], 
              label='Faulted', 
              linewidth='3')
    p[1].plot(FLTdiff[0], FLTdiff[1] + diffAdj, 
              color='#BEBEBE', 
              label='Difference', 
              linewidth='1')

    # set x-axis label
    x_label = r'Q (\AA" r"$^{-1}$, $\lambda=$' + str(wl) + r' \AA)'
    fig.supxlabel(x_label, fontsize=18)
    
    # set y-axis label
    if normalized == True:
        y_label = 'Intensity (counts, normalized)'
    elif normalized == False:
        y_label = 'Intensity (counts)'
    p[0].set_ylabel(y_label, fontsize=18)
    
    # legend labels
    UFlabel = '\n'.join(('Unfaulted', re.sub('x', str(nStacks), r'$N = x$')))
    FLTlabel = '\n'.join(('Faulted', re.sub('x', str(nStacks), r'$N = x$')))
    
    # legend handles
    obsHandle = mlines.Line2D([], [], label='Observed', 
                              color='white', marker='.', mfc='black', ms=15)
    UFhandle = mlines.Line2D([], [], label=UFlabel, color=simColors[0])
    FLThandle = mlines.Line2D([], [], label=FLTlabel, color=simColors[1])
    diffHandle = mlines.Line2D([], [], color='#BEBEBE', label='Difference')
    
    # add legend
    p[1].legend(handles=[obsHandle, UFhandle, FLThandle, diffHandle], 
                handlelength=1, fontsize='16', loc='upper left', 
                bbox_to_anchor=(1 + legendAdj[0], 1 + legendAdj[1]))
    
    # add fault parameters text box
    probText = re.sub('x', str(int(prob*100)), r'$P = x \%$')
    sVecText = r'$\vec{S} = \left[ x, y, z \right]$'
    varList = ['x', 'y', 'z']
    for i in range(3):
        subText = re.sub(varList[i], sVec[i], sVecText)
        sVecText = subText 
    params = '\n'.join((sVecText, probText))
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.5, pad=0.3)
    p[1].text(xLim[0] + boxAdj[0], yLim[1] + boxAdj[1], 
              params, bbox=props, ha='left', va='top', 
              fontsize='16', linespacing=2)
    
    plt.subplots_adjust(wspace=0.05)
    return(p)

#---------------------------------------------------------------------------------
# add (hkl) labels to plot -------------------------------------------------------
#---------------------------------------------------------------------------------
def addPeakLabels(plot, hkl, xPos, yPos, color=None, size=None):
    '''
    Parameters
    ----------
    plot
        Figure : variable name of figure
    hkl
        list (str) : text labels for (hkl) reflections of interest
    xPos
        list (float) : x-axis position for each text label
    yPos
        list (float) : y-axis position for each text label
    color
        str [optional] : label text color
        The default is black
    size
        str [optional] : label text size
        The default is 14
    '''
    if color is None:
        color='black'
    if size is None:
        size='14'
    
    for i in range(len(hkl)):
        plot.text(xPos[i], yPos[i], hkl[i], color=color, ha='center', 
                  va='center', fontsize=size)

#---------------------------------------------------------------------------------
# compare goodness of fit between data sets using difference of difference -------
# rows correspond to different models --------------------------------------------
# columns correspond to difference reflections of interest -----------------------
#---------------------------------------------------------------------------------
def fitCompare(rows, cols, diffQ, diffInts, xLims, yLim, wl, 
               rowLabels, colLabels, rowLabelAdj=None, colLabelAdj=None, 
               xLabelAdj=None, yLabelAdj=None, normalized=False, gradient=None):
    '''
    Parameters
    ----------
    rows
        int : number of rows
    cols 
        int : number of columns
    diffQ
        nparray : Q values (A^-1)
    diffInts
        list (nparray) : difference curve intensity values
        Format: list of lists where each row corresponds to list [c1, c2, ...]
    xLim
        list (float) : x-axis minimum and maximum, tuple
    yLim
        list (float) : y-axis minimum and maximum, tuple
    wl
        float : instrument wavelength (A)
    rowLabels
        list (str) : text labels for each row
    colLabels
        list (str) : text labels for each column
    rowLabelAdj
        float [optional] : adjustment parameter for row text position
        The default is 0.01
    colLabelAdj
        float [optional] : adjustment parameter for column text position
        The default is 0.01
    xLabelAdj
        float [optional] : adjustment parameter for x-axis label position
        The default is 0
    yLabelAdj
        float [optional] : adjustment parameter for y-axis label position
        The default is 0
    normalized
        bool [optional] : set to True if intensity values are normalized
        The default is False
    gradient
        list (str) [optional] : hex codes for corners of gradient map
        The default is ['#00C6BF', '#009AE1', '#5D7AD3', '#B430C2']
    '''
    
    fig, (p) = plt.subplots(rows, cols, figsize=(rows*2, cols))
    
    if rowLabelAdj is None:
        rowLabelAdj = 0.01
    if colLabelAdj is None:
        colLabelAdj = 0.01
    if xLabelAdj is None:
        xLabelAdj = 0
    if yLabelAdj is None:
        yLabelAdj = 0
    if gradient is None:
        gradient = ['#00C6BF', '#009AE1', '#5D7AD3', '#B430C2']
    
    g = gradientGen2D(gradient, rows, cols)
    
    # plot fit differences
    for row in range(rows):
        for col in range(cols):
            p[row][col].plot(diffQ, diffInts[row][col],
                             color = g[row][col].hex)
            # set axis limits
            p[row][col].set_xlim(xLims[col][0], xLims[col][1])
            p[row][col].set_ylim(yLim)
            p[row][col].tick_params(axis='both', labelsize='14')
  
    # format axes
    for row in range(rows-1):
        for col in range(cols):
            p[row][col].get_xaxis().set_visible(False)
    for row in range(rows):
        for col in range(1, cols):
            p[row][col].get_yaxis().set_visible(False)
    
    # set x-axis label
    x_label = r'Q (\AA' r'$^{-1}$, $\lambda=$' + str(wl) + r' \AA)'
    fig.supxlabel(x_label, fontsize=16, y=xLabelAdj)
    
    # set y-axis label
    if normalized == True:
        y_label = r'Diff$_{\mathrm{UF}} -$ Diff$_{\mathrm{F}}$ (counts, normalized)'
    elif normalized == False:
        y_label = r'Diff$_{\mathrm{UF}} -$ Diff$_{\mathrm{F}}$ (counts)'
    fig.supylabel(y_label, fontsize=16, x=yLabelAdj)
    
    # set plot labels
    y_mid = ((yLim[1] - yLim[0]) / 2) + yLim[0]
    x_end = xLims[-1][1]
    
    for row in range(rows):
        p[row][-1].text(x_end + rowLabelAdj, y_mid, 
                        rowLabels[row],
                        color=g[row][-1].hex, 
                        fontsize='16', ha='left', va='center')
    for col in range(cols):
        x_mid = ((xLims[col][1] - xLims[col][0]) / 2) + xLims[col][0]
        p[0][col].text(x_mid, yLim[1] + colLabelAdj, 
                       colLabels[col], 
                       color=g[0][col].hex, 
                       fontsize='16', ha='center', va='bottom')
        
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    return(p)

#---------------------------------------------------------------------------------
# set spacing of tick marks ------------------------------------------------------
#---------------------------------------------------------------------------------
def setTickSpacing(plot, axis, interval):
    '''
    Parameters
    ----------
    plot
        Figure : variable name of figure
    axis
        str : axis to adjust; set to 'x', 'y', or 'both'
    interval
        float : spacing interval
    '''
    if axis == 'x':
        plot.xaxis.set_major_locator(ticker.MultipleLocator(interval))
    elif axis == 'y':
        plot.yaxis.set_major_locator(ticker.MultipleLocator(interval))
    elif axis == 'both':
        plot.xaxis.set_major_locator(ticker.MultipleLocator(interval))
        plot.yaxis.set_major_locator(ticker.MultipleLocator(interval))


