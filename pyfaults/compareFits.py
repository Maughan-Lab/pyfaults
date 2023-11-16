##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import pandas as pd
import glob
import numpy as np

#---------------------------------------------------------------------------------
# generates a CSV file with information about broadened reflections --------------
#---------------------------------------------------------------------------------
def peakParams(labels, Qvals, savePath, pw=None):
    '''
    Parameters
    ----------
    labels
        list (str) : (hkl) labels for each peak
    Qvals
        list (float) : Q (A^-1) positions of each peak
    savePath
        str : file directory to export peak info as CSV
    pw
        float [optional] : broadened peak width (A^-1), default is 0.01
    '''
    
    peaks = pd.DataFrame()
    if pw is None:
        pw = 0.01
    
    Q_min = []
    Q_max = []
    # determine minimum and maximum Q for each peak
    for i in range(len(labels)):
        qMin = float('%.3f'%(Qvals[i] - (pw/2)))
        qMax = float('%.3f'%(Qvals[i] + (pw/2)))
        Q_min.append(qMin)
        Q_max.append(qMax)
        
    peaks['(hkl)'] = labels
    peaks['Q'] = Qvals
    peaks['Q min'] = Q_min
    peaks['Q max'] = Q_max
    
    peaks.to_csv(savePath + 'peak_info.csv', index=False)
    
    return

#---------------------------------------------------------------------------------
# calculates area under difference of fit residuals curve to compare models ------
#---------------------------------------------------------------------------------
def compareFits(fitDiffPath, peaksPath, peaksName):
    '''
    Parameters
    ----------
    fitDiffPath
        str : fit difference curve data file directory
    peaksPath
        str : peak info file directory
    peaksName
        str : file name of peak info CSV generated from peakParams
        
    Returns
    -------
    compareFits
        DataFrame : tabulated data of fit comparisons, includes the following as columns
            'Model' -- supercell model
            '(hkl)' -- broadened reflection
            'Peak Index (Q)' -- Q (A^-1) index for reflection
            'Faulted Model Fit' -- fit comparison of faulted model to unfaulted model, 
                either 'improved', 'worsened', or 'no change'
    '''
    import pyfaults as pf
    
    peaks = pd.read_csv(peaksPath + peaksName + '.csv')
    
    fitDiffs = glob.glob(fitDiffPath + '/*.txt')
    
    modelNames = []
    peakNames = []
    peakQvals = []
    fitResults = []
    
    # loop through each model
    for f in fitDiffs:
        removeExt = f.split('.')
        fn = removeExt[0].split('\\')
        saveName = fn[-1].split('_fitDiff')
        
        fitDiffQ, fitDiff = pf.importFile(fitDiffPath, fn[-1])
        
        # loop through each (hkl) index
        for p in peaks.index:
            qmin = peaks['Q min'][p]
            qmax = peaks['Q max'][p]
            
            # truncate fit diff curve
            truncFitDiff = []
            for q in range(len(fitDiffQ)):
                if fitDiffQ[q] >= qmin and fitDiffQ[q] <= qmax:
                    truncFitDiff.append(fitDiff[q])
                    
            # calculate area under curve
            diffArea = np.cumsum(truncFitDiff)
            if diffArea[-1] < 0:
                result = 'Improved'
            elif diffArea[-1] > 0:
                result = 'Worsened'
            elif diffArea[-1] == 0:
                result = 'No Change'
                
            modelNames.append(saveName[0])
            peakNames.append(peaks['(hkl)'][p])
            peakQvals.append(peaks['Q'][p])
            fitResults.append(result)
    
    # store results in data frame
    compareFits = pd.DataFrame()
    compareFits['Model'] = modelNames
    compareFits['(hkl)'] = peakNames
    compareFits['Peak Index (Q)'] = peakQvals
    compareFits['Faulted Model Fit'] = fitResults
    
    compareFits.to_csv(fitDiffPath + 'compareFitResults.csv', index=False)
    
    return compareFits

