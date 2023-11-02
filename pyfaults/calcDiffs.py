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
    '''
    import os, glob
    
    # create 'diffCurves' folder in file directory
    if os.path.exists(simPath + 'diffCurves/') == False:
        os.mkdir(simPath + 'diffCurves/')
    
    # import experimental data
    exptFile = exptName.split('.')
    ext = '.' + exptFile[1]
    exptQ, exptInts = pf.importExpt(exptPath, exptFile[0], wl, maxTT, ext=ext)
    
    # import simulation data
    sims = glob.glob(simPath + '/*.txt')
    
    for f in sims:
        removeExt = f.split('.')
        fn = removeExt[0].split('\\')
        
        # calculate expt vs model difference
        q, ints = pf.importFile(simPath, fn[-1])
        diffQ, diffInts = pf.diffCurve.diffCurve(exptQ, q, exptInts, ints)
        
        # save difference curve
        with open(simPath + 'diffCurves/' + fn[-1] + '_exptDiff.txt', 'w') as x:
            for (q, ints) in zip(diffQ, diffInts):
                x.write('{0} {1}\n'.format(q, ints))
            x.close()
    
    return