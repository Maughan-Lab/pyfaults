########################################################################
# Author: Sinclair R. Combs
########################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
rc('text', usetex=True)
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']},size='14')
rc('text.latex',preamble=r'\usepackage{sfmath}')

''' AUTOSEARCH METHODS '''
#-----------------------------------------------------------------------
# generate base unit cell ----------------------------------------------
#-----------------------------------------------------------------------
def setUnitCell(csvPath, csvFile, lattice, layerList, stackDir,
                childLayerList=None, childGenVecList=None):
    '''
    Parameters
    ----------
    see pyfaults.unitcell.genUnitCell()
    '''
    import pyfaults.unitcell as ucell
    return ucell.genUnitCell("Base Stack, N=1", csvPath, csvFile, 
                      lattice, layerList, stackDir, 
                      childLayerList, childGenVecList)

#-----------------------------------------------------------------------
# generate data frame with broadened peak info -------------------------
#-----------------------------------------------------------------------
def peakParams(peakLabels, peakQ, calcPW):
    '''
    Parameters
    ----------
    peakLabels
        list (str) : (hkl) labels for each peak
    peakQ
        list (float) : Q (A^-1) positions of each peak
    calcPW
        float : broadened peak width (A^-1)

    Returns
    -------
    peakDF
        DataFrame : tabulated peak data, includes the following as columns
            'Reflection' -- (hkl) labels
            'Q' -- Q indexes (A^-1)
            'Peak Range' -- minimum and maximum Q values for each peak (A^-1)
    '''
    peakDF = pd.DataFrame()
    peakRanges = []
    # determine minimum and maximum Q for each peak
    for i in range(len(peakLabels)):
        qMin = float('%.3f'%(peakQ[i] - (calcPW/2)))
        qMax = float('%.3f'%(peakQ[i] + (calcPW/2)))
        qRange = [qMin, qMax]
        peakRanges.append(qRange)
    # add entry for each peak
    peakDF['Reflection'] = peakLabels
    peakDF['Q'] = peakQ
    peakDF['Peak Range'] = peakRanges
    return peakDF

#-----------------------------------------------------------------------
# compare fits to broadened peaks between models -----------------------
#-----------------------------------------------------------------------
def compareFits(Q, fitDiffDF, peakDF):
    '''
    Parameters
    ----------
    Q
        nparray : Q-values (A^-1)
    fitDiffDF
        DataFrame : tabulated data, includes the following
            'Model' -- unique identifier
            'Stacking Vector' -- stacking vector in string format
            'Stacking Probability' -- stacking fault probability in str format
            'S_x' -- x-component of stacking vector
            'S_y' -- y-component of stacking vector
            'S_z' -- z-component of stacking vector
            'P' -- fault probability (0 to 1)
            'Simulated Q' -- calculated Q values (A^-1)
            'Simulated Intensity' -- calculated intensity values
            'Expt vs. Model Q' -- difference curve Q Values (A^-1)
            'Expt vs. Model Difference' -- difference curve intensity values
            'UF vs. FLT Model' -- intensity values of fit difference curve
    peakDF
        DataFrame : broadened peak info

    Returns
    -------
    fitComp
        list (str, DataFrame) : each list entry contains the model name and a
        dataframe including the following
            'Peak' -- (hkl) identifier
            'UF vs. FLT Fit Difference' -- fit difference intensity values
            'Sum' -- area under the fit difference curve
            'Fit Result' -- indicator of change in fit between unfaulted and
            faulted models ('improved', 'worsened', or 'unchanged')
    '''
    
    fitComp = []
    # for loop over each supercell
    for i in fitDiffDF.index:
        name = fitDiffDF['Model'][i]
        compDF = pd.DataFrame()
        Q = fitDiffDF['Expt vs. Model Q'][i]
        fullFitDiff = fitDiffDF['UF vs. FLT Model'][i]
        
        truncFitCol = []
        sumCol = []
        resultCol = []
        # unfaulted supercell entry
        if i == 0:
            truncFitCol.append(np.nan)
            sumCol.append(np.nan)
            resultCol.append(np.nan)
        # faulted supercell entries
        if i > 0:
            Q = fitDiffDF['Expt vs. Model Q'][i]
            fullFitDiff = fitDiffDF['UF vs. FLT Model'][i]
            
            # for loop over each broadened reflection
            for j in peakDF.index:
                qRange = peakDF['Peak Range'][j]
                qMin = qRange[0]
                qMax = qRange[1]
                
                # truncate according to peak range
                truncFitDiff = []
                for val in range(len(Q)):
                    if Q[val] >= qMin and Q[val] <= qMax:
                        truncFitDiff.append(fullFitDiff[val])
                        
                # determine area under curve and assess change in fit
                diffSum = np.cumsum(truncFitDiff)
                if diffSum[-1] < 0:
                    result = 'Improved'
                elif diffSum[-1] > 0:
                    result = 'Worsened'
                elif diffSum[-1] == 0:
                    result = 'No Change'
                truncFitCol.append(truncFitDiff)
                sumCol.append(diffSum)
                resultCol.append(result)
            
            # add entry for each broadened reflection
            compDF['Peak'] = peakDF['Reflection']
            compDF['UF vs. FLT Fit Difference'] = truncFitCol
            compDF['Sum'] = sumCol
            compDF['Fit Result'] = resultCol
            # create list entry of model name with results
            fitComp.append([name, compDF])
        
        return fitComp

#-----------------------------------------------------------------------
# automated parameter grid searching -----------------------------------
#-----------------------------------------------------------------------
def autoSearch(savePath, unitcell, expt, nStacks, fltLayer, 
               probList, sVecList, peakDF, wl, maxTT, simPW=None):
    '''
    Parameters
    ----------
    savePath
        str : file save directory
    unitcell
        Unitcell : base unit cell stack
    expt
        list (nparray) : experimental data as [Q, ints]
    nStacks
        int : number of stacks per supercell
    fltLayer
        str : name of faulted layer
    probList
        list (float) : stacking fault probabilities to search
    sVecList
        list (nparray) : stacking vectors to search
    peakDF
        DataFrame : broadened peak info
    wl
        float : instrument wavelength (A)
    maxTT
        float : maximum 2theta value (degrees)
    simPW
        float [optional]: simulated peak broadening term (applied to all peaks)
        The default is 0.001

    Returns
    -------
    fitComp
        list (str, DataFrame) : each list entry contains the model name and a
        dataframe including the following
            'Peak' -- (hkl) identifier
            'UF vs. FLT Fit Difference' -- fit difference intensity values
            'Sum' -- area under the fit difference curve
            'Fit Result' -- indicator of change in fit between unfaulted and
            faulted models ('improved', 'worsened', or 'unchanged')
    fitDiffDF
        DataFrame : tabulated data, includes the following
            'Model' -- unique identifier
            'Stacking Vector' -- stacking vector in string format
            'Stacking Probability' -- stacking fault probability in str format
            'S_x' -- x-component of stacking vector
            'S_y' -- y-component of stacking vector
            'S_z' -- z-component of stacking vector
            'P' -- fault probability (0 to 1)
            'Simulated Q' -- calculated Q values (A^-1)
            'Simulated Intensity' -- calculated intensity values
            'Expt vs. Model Q' -- difference curve Q Values (A^-1)
            'Expt vs. Model Difference' -- difference curve intensity values
            'UF vs. FLT Model' -- intensity values of fit difference curve
    '''
    import pyfaults.supercell as sup
    import pyfaults.simXRD as sim
    import pyfaults.diffCalc as diff
    
    # generate supercells ----------------------------------------------
    cellDF = sup.genSupercells(unitcell, nStacks, fltLayer, 
                               probList, sVecList, savePath)
    print('Finished generating supercells')
    print('See ' + savePath + 'CIFs')
    
    # simulate XRD -----------------------------------------------------
    if simPW is None:
        simPW = 0.001
    simDF = sim.calcSims(cellDF, savePath, wl, maxTT, simPW)
    print('Finished simulating X-ray diffraction patterns')
    print('See ' + savePath + 'sims')
    
    # calculate difference curves --------------------------------------
    exptDiffDF = diff.calcDiffs(savePath, simDF, expt)
    print('Finished calculating difference curves')
    print('See ' + savePath + 'diffCurves')
    
    # calculate fit differences ----------------------------------------
    fitDiffDF = diff.calcFitDiffs(savePath, exptDiffDF)
    print('Finished calculating fit difference curves')
    print('See ' + savePath + 'fitDiffCurves')
    
    # compare peak fits ------------------------------------------------
    # set Q values
    Q = fitDiffDF['Simulated Q'][0]
    
    fitComp = compareFits(Q, fitDiffDF, peakDF)
    
    # export fit comparisons
    with open (savePath, 'w') as f:
        for i in range(len(fitComp)):
            modelResult = fitComp[i][1]
            f.write('Model: ' + fitComp[i][0] + '\n')
            for j in modelResult.index:
                f.write('Peak: ' + modelResult['Peak'][j] + ' ' +\
                        modelResult['Fit Result'][j] + '\n')
            f.write('\n')

    return fitComp, fitDiffDF

#-----------------------------------------------------------------------
# plots peakFitcompare() results ---------------------------------------
#-----------------------------------------------------------------------
def plotSearchResults(fitDiffDF, peakDF, xSpacing, wl, yLim, 
                      rowAdj=0.0, colAdj=0.0, 
                      xLabelAdj=0.0, yLabelAdj=0.0):
    '''
    Parameters
    ----------
    fitDiffDF
        DataFrame : tabulated data returned from calcFitDiffs()
    peakDF
        DataFrame : tabulated data returned from makePeakDF()
    xSpacing
        float : adjustment parameter for spacing around peaks (x-axis units)
    wl
        float : instrument wavelength (A)
    yLim
        list (flaot) : y-axis minimum and maximum
    rowAdj
        float [optional] : adjustment parameter for row label x-axis positions
        The default is 0.0
    colAdj
        float [optional] : adjustment parameter for col label y-axis positions
        The default is 0.0
    xLabelAdj
        float [optional] : adjustment parameter for x-axis label position
        The default is 0.0
    yLabelAdj
        float [optional] : adjustment parameter for y-axis label position
        The default is 0.0

    Returns
    -------
    p : Figure
    '''
    from pyfaults.plotXRD import gradientGen2D
    
    # set number of rows and columns
    rows = len(fitDiffDF)-1
    cols = len(peakDF)
    
    # generate color map
    g = gradientGen2D(['#00C6BF', '#009AE1', '#5D7AD3', '#B430C2'], rows, cols)
    
    # create new Figure instance
    fig, (p) = plt.subplots(rows, cols, figsize=(rows*2, cols))
    # set Q values
    Q = fitDiffDF['Expt vs. Model Q'][0]
    
    # plot fit difference curves
    for r in range(rows):
        for c in range(cols):
            ints = fitDiffDF['UF vs. FLT Model'][r+1]
            p[r][c].plot(Q, ints, color=g[r][c].hex)
            
    # set x-axis limits for each column (broadened peak)    
    xLims = []
    for i in peakDF.index:
        qRange = peakDF['Peak Range'][i]
        qMin = qRange[0] - xSpacing
        qMax = qRange[1] + xSpacing
        xLims.append([qMin, qMax])
    # set y-axis limits
    for r in range(rows):
        for c in range(cols):
            p[r][c].set_xlim(xLims[c][0], xLims[c][1])
            p[r][c].set_ylim(yLim)
            p[r][c].tick_params(axis='both', labelsize='14')
            if r != (rows-1):
                p[r][c].get_xaxis().set_visible(False)
            if c != 0:
                p[r][c].get_yaxis().set_visible(False)
    
    # add axis labels
    xLabel = r'Q (\AA' r'$^{-1}$, $\lambda=$' + str(wl) + r' \AA)'
    fig.supxlabel(xLabel, fontsize=16, y=xLabelAdj)
    yLabel = r'Diff$_{\mathrm{UF}} -$ Diff$_{\mathrm{F}}$ (counts, normalized)'
    fig.supylabel(yLabel, fontsize=16, x=yLabelAdj)
    
    # set row and column labels
    rowLabels = []
    for r in range(rows):
        sVec = str(fitDiffDF['Stacking Vector'][r+1])
        prob = str(fitDiffDF['Stacking Probability'][r+1])
        label = '\n'.join((sVec, prob))
        rowLabels.append(label)
    colLabels = peakDF['Reflection']
    
    # add row and column labels
    yMid = ((yLim[1] - yLim[0]) / 2) + yLim[0]
    xEnd = xLims[-1][1]
    for r in range(rows):
        p[r][-1].text(xEnd+rowAdj, yMid, rowLabels[r], color=g[r][-1].hex, 
                      fontsize='14', ha='left', va='center')
    for c in range(cols):
        xMid = ((xLims[c][1] - xLims[c][0]) / 2) + xLims[c][0]
        p[0][c].text(xMid, yLim[1]+colAdj, colLabels[c], color=g[0][c].hex, 
                     fontsize='14', ha='center', va='bottom')
    
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    return(p)
















        
        