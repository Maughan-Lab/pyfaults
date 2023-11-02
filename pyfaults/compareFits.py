########################################################################
# Author: Sinclair R. Combs
########################################################################
import pandas as pd

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

    Returns
    -------
    peaks
        DataFrame : tabulated peak data, includes the following as columns
            '(hkl)' -- (hkl) plane corresponding to reflection
            'Q' -- Q indexes (A^-1)
            'Q min' -- lower Q limit of broadened peak
            'Q max' -- upper Q limit of broadened peak
    '''
    peaks = pd.DataFrame()
    if pw is None:
        pw = 0.01
    
    Q_min = []
    Q_max = []
    # determine minimum and maximum Q for each peak
    for i in range(len(labels)):
        qMin = float('%.3f'%(peakQ[i] - (pw/2)))
        qMax = float('%.3f'%(peakQ[i] + (pw/2)))
        Q_min.append(qMin)
        Q_max.append(qMax)
        
    peaks['(hkl)'] = labels
    peaks['Q'] = Qvals
    peaks['Q min'] = Q_min
    peaks['Q max'] = Q_max
    
    peaks.to_csv(savePath + 'peak_info.csv')
    
    return peaks

def compareFits(fitDiffPath, fitDiffCSV, peaksPath, peaksName):
    '''
    Parameters
    ----------
    fitDiffPath
        str : fit difference curve data file directory
    fitDiffCSV
        str : file name of metadata CSV generated from calcFitDiffs method
    peaksPath
        str : peak info file directory
    peaksName
        str : file name of peak info CSV generated from peakParams

    Returns
    -------
    compareFits
        DataFrame : tabulated fit comparison data, includes the following as columns
            'Model' -- unique identifier
            '(hkl)' -- (hkl) plane corresponding to reflection
            'Peak Index (Q)' -- Q indexes (A^-1)
            'Faulted Model Fit' -- result of fit comparison for each peak
    '''
    
    fitDiffDF = pd.read_csv(fitDiffPath + fitDiffCSV + '.csv')
    peaks = pd.read_csv(peaksPath + peaksName + '.csv')
    
    modelNames = []
    peakNames = []
    peakQvals = []
    fitResults = []
    
    # loop through each supercell
    for i in fitDiffDF.index:
        name = fitDiffDF['Model'][i]
        Q = fitDiffDF['Q'][i]
        fitDiff = fitDiffDF['UF vs FLT'][i]
        
        # skip unfaulted supercell
        if i > 0:
            # loop through each (hkl) index
            for j in peaks.index:
                qmin = peaks['Q min'][j]
                qmax = peaks['Q max'][j]
                
                # truncate fit diff curve
                truncFitDiff = []
                for q in range(len(Q)):
                    if Q[q] >= qmin and Q[q] <= qmax:
                        truncFitDiff.append(fitDiff[q])
                        
                # calculate area under curve
                diffArea = np.cumsum(truncFitDiff)
                if diffArea[-1] < 0:
                    result = 'Improved'
                elif diffArea[-1] > 0:
                    result = 'Worsened'
                elif diffArea[-1] == 0:
                    result = 'No Change'
                    
                modelNames.append(name)
                peakNames.append(peaks['(hkl)'][j])
                peakQvals.append(peaks['Q'][j])
                fitResults.append(result)
    
    # store results in data frame
    compareFits = pd.DataFrame()
    compareFits['Model'] = modelNames
    compareFits['(hkl)'] = peakNames
    compareFits['Peak Index (Q)'] = peakQvals
    compareFits['Faulted Model Fit'] = fitResults
    
    compareFits.to_csv(fitDiffPath + 'compareFitResults.csv')
    
    return compareFits

