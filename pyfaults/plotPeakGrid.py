##################################################################################
# Author: Sinclair R. Combs
##################################################################################

from colour import Color
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.pyplot import rc
rc('text', usetex=True)
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']},size='14')
rc('text.latex',preamble=r'\usepackage{sfmath}')


########################################################################################
# 2D GRADIENT GENERATOR
########################################################################################
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


########################################################################################
# PEAK GRID PLOTTING FUNCTION
########################################################################################
def plotPeakGrid(rows, cols, qVals, intsVals, colLims, expt=None, figSize=None, 
                 gradient=None, rowLabels=None, rowAdj=0, colLabels=None, colAdj=0):
    
    if figSize == None:
        figSize = (8,8)
    if gradient is None:
        gradient = ['#00C6BF', '#009AE1', '#5D7AD3', '#B430C2']
    g = gradientGen2D(gradient, rows, cols)
    
    fig, (p) = plt.subplots(rows, cols, figsize=(figSize))

    # plot data
    for row in range(rows):
        for col in range(cols):
            if expt is not None:
                exptMin = np.min(expt[1])
                p[row][col].scatter(expt[0], expt[1], color='black', marker='.', s=8)
            else:
                exptMin = 0
            p[row][col].plot(qVals[row], intsVals[row]+exptMin, color = g[row][col].hex)
            
    # determine y-axis limits
    yLims = (0, 0)
    for i in range(len(intsVals)):
        currMin = np.min(intsVals[i])
        currMax = np.max(intsVals[i])
        if currMin < yLims[0]:
            yLims[0] = currMin
        if currMax > yLims[1]:
            yLims[1] = currMax
            
    # set axis limits
    for row in range(rows):
        for col in range(cols):
            p[row][col].set_xlim(colLims[col][0], colLims[col][1])
            p[row][col].set_ylim(yLims)
            p[row][col].tick_params(axis='both', labelsize='14')
            
            if row < (rows-1):
                p[row][col].get_xaxis().set_visible(False)
            if col > 0:
                p[row][col].get_yaxis().set_visible(False)
                
    # set row labels
    yMid = ((yLims[1] - yLims[0]) / 2) + yLims[0]
    if rowLabels is not None:
        for row in range(rows):
            p[row][-1].text(colLims[-1][1] + rowAdj, yMid, rowLabels[row], 
                            color='black', fontsize='16', ha='left', va='center')
            
    # set column labels
    if colLabels is not None:
        for col in range(cols):
            xMid = ((colLims[col][1]-colLims[col][0]) / 2) + colLims[col][0]
            p[0][col].text(xMid, yLims[1] + colAdj, colLabels[col], 
                           color='black', fontsize='16', ha='left', va='center')
    
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    return(p)

    
########################################################################################
# PLOT EDITING FUNCTIONS
########################################################################################

#---------------------------------------------------------------------------------------
# set spacing of tick marks ------------------------------------------------------------
#---------------------------------------------------------------------------------------
def setTickSpacing(p, axis, interval):
    '''
    Parameters
    ----------
    p
        Figure : plot to edit
    axis
        str : axis to adjust; set to 'x', 'y', or 'both'
    interval
        float : spacing interval
    '''
    if axis == 'x':
        p.xaxis.set_major_locator(ticker.MultipleLocator(interval))
    elif axis == 'y':
        p.yaxis.set_major_locator(ticker.MultipleLocator(interval))
    elif axis == 'both':
        p.xaxis.set_major_locator(ticker.MultipleLocator(interval))
        p.yaxis.set_major_locator(ticker.MultipleLocator(interval))

#---------------------------------------------------------------------------------------
# set x-axis label ---------------------------------------------------------------------
#---------------------------------------------------------------------------------------
def setXLabel(p, xLabel, wl, size=None):
    '''
    Parameters
    ----------
    p
        Figure : plot to edit
    xLabel
        str : x-axis label, set to 'Q' to use default Q label
    wl
        float : instrument wavelength (A)
    size
        int (optional): font size, default is 16
    '''
    if size == None:
        size = '16'

    if xLabel == 'Q':
        txt = r'$Q$ (\AA $^{-1}$, $\lambda=$' + str(wl) + r'\AA)'
        if p == 'fig':
            p.supxlabel(txt, fontsize=size)
        else:
            p.set_xlabel(txt, fontsize=size)
    else:
        p.set_xlabel(xLabel, fontsize=size)
    return

#---------------------------------------------------------------------------------------
# set y-axis label ---------------------------------------------------------------------
#---------------------------------------------------------------------------------------
def setYLabel(p, yLabel, labelType, size=None):
    '''
    Parameters
    ----------
    p
        Figure : plot to edit
    yLabel
        str : y-axis label
    labelType
        str : 'ints', 'diff', or 'fitDiff'
    size
        int (optional): font size, default is 16
    '''
    if size == None:
        size = '16'
        
    if labelType == 'ints':
        txt  = 'Intensity (a.u.)'
    elif labelType == 'diff':
        txt = r'$\mathrm{D}_{\mathrm{obs-sim}} = \mathrm{I}_{\mathrm{obs}} - \mathrm{I}_{\mathrm{sim}}$ (a.u.)'
    elif labelType == 'fitDiff':
        txt = r'$mathrm{R}_\mathrm{D} = \mathrm{D}_{\mathrm{obs-ideal}} - \mathrm{D}_{\mathrm{obs-fault}}$ (a.u.)'
        
    if p == 'fig':
        p.supylabel(txt, fontsize=size)
    else:
        p.set_ylabel(txt, fontsize=size)
    return

#---------------------------------------------------------------------------------------
# change axis label size ---------------------------------------------------------------
#---------------------------------------------------------------------------------------           
def adjAxSize(p, size):
    '''
    Parameters
    ----------
    p
        Figure : plot to edit
    size
        int : font size, default is 16
    '''
    p.tick_params(axis='both', labelsize=size)
    return