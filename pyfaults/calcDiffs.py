##################################################################################
# Author: Sinclair R. Combs
##################################################################################

def calcDiffs(exptPath, exptName, simPath, wl, maxTT):
    '''
    Parameters
    ----------
    exptPath
        str : experimental XRD file directory
    exptName
        str : file name of experimental data, including extension
    simPath
        str : simulation data file directory
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
    import os, glob
    import pandas as pd
    
    # create 'diffCurves' folder in file directory
    if os.path.exists(simPath + 'diffCurves/') == False:
        os.mkdir(simPath + 'diffCurves/')
    
    # import experimental data
    exptFile = exptName.split('.')
    ext = '.' + exptFile[1]
    exptQ, exptInts = pf.importExpt(exptPath, exptFile[0], wl, maxTT, ext=ext)
    
    # import simulation data
    sims = glob.glob(simPath + '/.txt')
    
    for f in sims:
        # calculate expt vs model difference
        q, ints = pf.importFile(simPath, f)
        diffQ, diffInts = pf.diffcurve.diffCurve(exptQ, q, exptInts, ints)
        
        # save difference curve
        with open(simPath + 'diffCurves/' + f + '_exptDiff.txt', 'w') as x:
            for (q, ints) in zip(diffQ, diffInts):
                x.write('{0} {1}\n'.format(q, ints))
            x.close()
    
    return diffData