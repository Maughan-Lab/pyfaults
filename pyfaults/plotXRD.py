########################################################################################
# Author: Sinclair R. Combs
########################################################################################

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


########################################################################################
# PLOTTING FUNCTIONS
########################################################################################

#---------------------------------------------------------------------------------------
# generates single plot ----------------------------------------------------------------
#---------------------------------------------------------------------------------------
def singlePlot(qVals, intsVals, expt=None, color=None, label=None):
    '''
    Parameters
    ----------
    qVals
        nparray : list of Q values (A^1)
    intsVals
        nparray : list of intensity values
    expt
        list (nparray) (optional) : experimental Q and intensity values
    color
        str (optional) : plot line color, default is '#00C6BF'
    label
        str (optional) : plot line label for legend, no legend if set to None
    '''
    
    fig, (p) = plt.subplots(1, figsize=(8, 8))
    
    # set color
    if color is None:
        color = '#00C6BF'
    
    # plot experimental data
    if expt is not None:
        exptMin = np.min(expt[1])
        p.scatter(expt[0], expt[1], color='black', label='Observed', marker='.', s=8)
    else:
        exptMin = 0
        
    # plot data
    p.plot(qVals, intsVals+exptMin, color=color, linewidth='2.5', label=label)
    
    # add legend
    if label is not None:
        p.legend(handlelength=1, fontsize='16')
    
    return(p)

#---------------------------------------------------------------------------------------
# generates stacked plot ---------------------------------------------------------------
#---------------------------------------------------------------------------------------
def stackedPlot(qVals, intsVals, expt=None, labels=None, gradient=None, labelPos=(0,0)):
    '''
    Parameters
    ----------
    qVals
        list (nparray) : list of Q value datasets
    intsVals
        list (nparray) : list of intensity value datasets
    expt
        list (nparray) (optional) : experimental Q and intensity values
    gradient
        list (str) (optional) : first and last gradient color hex codes
        The default is ['#00C6BF', '#B430C2']
    labels
        list (str) (optional) : list of text labels for each dataset
        The default is None
    labelPos
        list (float) (optional) : text label position
        The default is (0,0)
    '''
    num = len(qVals)
    
    fig, (p) = plt.subplots(num, 1, figsize=(8, 8))
    
    if gradient is None:
        gradient = ['#00C6BF', '#B430C2']
    
    g = list(Color(gradient[0]).range_to(Color(gradient[1]), num))
    
    # plot experimental data
    if expt is not None:
        exptMin = np.min(expt[1])
        for i in range(num):
            p[i].scatter(expt[0], expt[1], color='black', label='Observed', marker='.', s=8)
    else:
        exptMin = 0
        
    # plot stacked data
    for i in range(num):
        p[i].plot(qVals[i], intsVals[i]+exptMin, color=g[i].hex, linewidth='2')
    
    # add additional labels 
    if labels is not None:
        for i in range(num):
            p[i].text(labelPos[0], labelPos[1], labels[i], color=g[i].hex, fontsize='16', ha='right', va='top')
    
    # add legend
    p[0].legend(handlelength=1, fontsize='14')
    
    plt.subplots_adjust(hspace=0.05) 
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
# set axis limits ----------------------------------------------------------------------
#---------------------------------------------------------------------------------------
def setLimits(p, xLim=None, yLim=None, size=None, hide=False):
    '''
    Parameters
    ----------
    p
        Figure : plot to edit
    xLim
        float (optional) : x-axis minimum and maximum (min, max)
    yLim
        float (optional) : y-axis minimum and maximum (min, max)
    size
        int (optional): font size, default is 14
    hide
        bool (optional) : hide axis values, default is False
    '''
    if size == None:
        size = '14'

    if xLim is not None:
        p.set_xlim(xLim[0], xLim[1])
        p.tick_params(axis='x', labelsize=size)
        if hide == True:
            p.get_xaxis().set_visible(False)
    if yLim is not None:
        p.set_ylim(yLim[0], yLim[1])
        p.tick_params(axis='y', labelsize=size)
        if hide == True:
            p.get_yaxis().set_visible(False)
    return

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
def setYLabel(p, yLabel, size=None, norm=False):
    '''
    Parameters
    ----------
    p
        Figure : plot to edit
    yLabel
        str : y-axis label
    size
        int (optional): font size, default is 16
    norm
        bool (optional) : set to True if intensity is normalized, default is False
    '''
    if size == None:
        size = '16'
        
    if norm == False:
        txt = 'Intensity (a.u.)'
    elif norm == True:
        txt = 'Intensity (a.u., normalized)'
        
    if p == 'fig':
        p.supylabel(txt, fontsize=size)
    else:
        p.set_ylabel(txt, fontsize=size)
    return

#---------------------------------------------------------------------------------------
# add (hkl) labels ---------------------------------------------------------------------
#---------------------------------------------------------------------------------------
def addPeakLabels(p, peaks, size=None):
    '''
    Parameters
    ----------
    plot
        Figure : plot to edit
    peaks
        list (str, float): list of peak information, each entry should contain
        the (hkl) text label (str), the x-axis position (float), and the y-axis 
        position for the label (float) (i.e., ['(hkl)', x, y])
    size
        str (optional) : label text size, the default is 14
    '''
    if size is None:
        size='14'
    for i in range(len(peaks)):
        p.text(peaks[i][1], peaks[i][2], peaks[i][0], color='black', ha='center', va='center', fontsize=size)
    return

#---------------------------------------------------------------------------------------
# add difference curve -----------------------------------------------------------------
#---------------------------------------------------------------------------------------
def addDiffCurve(p, diff_q, diff_ints, offset=None, color=None):
    '''
    Parameters
    ----------
    p
        Figure : plot to edit
    diff_q
        nparray : difference curve Q values (A^-1)
    diff_ints
        nparray : difference curve intensity values
    offset
        float (optional) : vertical offset from data, default is -0.5
    color
        str (optional) : hex code of line color, default is '#BEBEBE'
    '''
    if offset == None:
        offset = -0.5
    if color == None:
        color = '#BEBEBE'
    
    p.plot(diff_q, diff_ints + offset, color=color, linewidth='2', label='Difference')
    return

#---------------------------------------------------------------------------------------
# add text box with faulting parameters (P, S) -----------------------------------------
#---------------------------------------------------------------------------------------       
def addFaultParams(p, prob, s, pos, color=None, size=None):
    '''
    Parameters
    ----------
    p
        Figure : plot to edit
    prob
        float : fault probability
    s
        nparray : stacking vector [x, y, z]
    pos
        nparray (float) (optional) : position of text box (x, y)
    color
        str (optional) : hex code of box fill color, default is white
    size
        str (optional) : text size, the default is 16
    '''
    pText = re.sub('x', str(int(prob*100)), r'$P = x \%$')
    sText = r'$\vec{S} = \left[ x, y, z \right]$'
    varList = ['x', 'y', 'z']
    for i in range(3):
        subText = re.sub(varList[i], s[i], sText)
        sText = subText 
    params = '\n'.join((sText, pText))
    
    if color == None:
        color = 'white'
    if size == None:
        size = '16'
    props = dict(boxstyle='square', facecolor=color, alpha=0.5, pad=0.3)
    
    p.text(pos[0], pos[1], params, bbox=props, ha='left', va='top', fontsize=size, linespacing=2)
    return

#---------------------------------------------------------------------------------------
# create legend handle -----------------------------------------------------------------
#---------------------------------------------------------------------------------------    
def createLegendHandle(label, color, handleType):
    '''
    Parameters
    ----------
    label
        str : handle text label
    color
        str : hex code of handle color
    handleType
        str : 'line' or 'dot'

    Returns
    -------
    handle
        Line2D : handle, can parse to function to add legend to plot
    '''
    if handleType == 'line':
        handle = mlines.Line2D([], [], label=label, color=color)
    elif handleType == 'dot':
        handle = mlines.Line2D([], [], label=label, color='white', marker='.', mfc=color, ms=15)
    
    return handle