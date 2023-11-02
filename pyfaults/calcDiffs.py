##################################################################################
# Author: Sinclair R. Combs
##################################################################################

def calcDiffs(exptPath, exptName, simPath, CSVName, wl, maxTT):
    '''
    Parameters
    ----------
    exptPath
        str : experimental XRD file directory
    exptName
        str : file name of experimental data, including extension
    simPath
        str : simulation data file directory
    CSVName
        str : file name of metadata CSV generated from calcSims method
    wl
        float : instrument wavelength (A)
    maxTT
        float : maximum 2theta value (degrees)

    Returns
    -------
    diffData
        DataFrame : tabulated diff curve data, includes the following as columns
            'Model' -- unique identifier
            'Expt vs. Model Diff Q' -- difference curve Q values
            'Expt vs. Model Diff Intensity' -- intensity differences between
            experimental XRD and simulated XRD
    '''
    import pyfaults as pf
    import os
    import pandas as pd
    
    # create 'diffCurves' folder in file directory
    if os.path.exists(simPath + 'diffCurves/') == False:
        os.mkdir(simPath + 'diffCurves/')
        
    exptFile = exptName.split('.')
    ext = '.' + exptFile[1]
    
    exptQ, exptInts = pf.importExpt(exptPath, exptFile[0], wl, maxTT, ext=ext)
    
    simDF = pd.read_csv(simPath + CSVName + '.csv')

    exptDiffQ = []
    exptDiffInts = []
    # calculate experimental vs model difference for each supercell
    for row in simDF.index:
        name = simDF['Model'][row]
        simQ = simDF['Sim Q'][row]
        simInts = simDF['Sim Norm Intensity'][row]
        diffQ, diffInts = pf.diffCurve.diffCurve(exptQ, simQ, exptInts, simInts)
        exptDiffQ.append(diffQ)
        exptDiffInts.append(diffInts)
        # save difference data
        with open(simPath + 'diffCurves/' + name + '_exptDiff.txt', 'w') as f:
            for (q, ints) in zip(diffQ, diffInts):
                f.write('{0} {1}\n'.format(q, ints))
        f.close()
        
    # store difference curve values in data frame
    diffData = pd.dataFrame()
    diffData['Model'] = simDF['Model']
    diffData['Expt vs. Model Diff Q'] = exptDiffQ
    diffData['Expt vs. Model Diff Intensity'] = exptDiffInts
    
    diffData.to_csv(simPath + 'diffCurves/diffCurve_info.csv')
    
    return diffData