##################################################################################
# Author: Sinclair R. Combs
##################################################################################

def calcFitDiffs(simPath, diffPath, CSVName):
    '''
    Parameters
    ----------
    simPath
        str : simulation data file directory
    diffPath
        str : difference curve data file directory
    CSVName
        str : file name of metadata CSV generated from calcDiffs method

    Returns
    -------
    fitDiffData
        DataFrame : tabulated fit diff data, includes the following as columns
            'Model' -- unique identifier
            'Q' -- Q values
            'UF vs FLT' -- intensity differences between unfaulted model fit and
            faulted model fit to experimental XRD
    '''
    import os
    import pandas as pd
    import numpy as np
    
    # create 'fitDiffCurves' folder in file directory
    if os.path.exists(simPath + 'fitDiffCurves/') == False:
        os.mkdir(simPath + 'fitDiffCurves/')
        
    diffDF = pd.read_csv(diffPath + CSVName + '.csv')
    
    Q = diffDF['Expt vs. Model Diff Q'][0]
    
    UFdiff = diffDF['Expt vs. Model Diff Intensity'][0]
    
    fitQ = []
    fitDiffs = []
    # calculate fit differences
    for i in diffDF.index:
        if i == 0:
            fitDiffs.append('None')
        elif i > 0:
            name = diffDF['Model'][i]
            modelDiff = diffDF['Expt vs. Model Diff Intensity'][i]
            diff = np.subtract(UFdiff, modelDiff)
            fitDiffs.append(diff)
            fitQ.append(Q)
            # save fit difference data
            with open(simPath + 'fitDiffCurves/' + name + '_fitDiff.txt', 'w') as f:
                for q, diff in range(len(Q, fitDiffs)):
                    f.write('{0} {1}\n'.format(q, diff))
            f.close() 
            
    # store fit differences in data frame
    fitDiffData = pd.dataFrame()
    fitDiffData['Model'] = diffDF['Model']
    fitDiffData['Q'] = fitQ
    fitDiffData['UF vs FLT'] = fitDiffs
    
    fitDiffData.to_csv(simPath + 'fitDiffCurves/fitDiff_info.csv')
    
    return fitDiffData