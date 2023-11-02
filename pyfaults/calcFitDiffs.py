##################################################################################
# Author: Sinclair R. Combs
##################################################################################

def calcFitDiffs(simPath, diffPath):
    '''
    Parameters
    ----------
    simPath
        str : simulation data file directory
    diffPath
        str : difference curve data file directory
    '''
    import pyfaults as pf
    import os, glob
    import numpy as np
    
    # create 'fitDiffCurves' folder in file directory
    if os.path.exists(simPath + 'fitDiffCurves/') == False:
        os.mkdir(simPath + 'fitDiffCurves/')
    
    # import difference curves
    unfaultedDiffQ, unfaultedDiff = pf.importFile(diffPath, 'Unfaulted_sim_exptDiff')
    diffCurves = glob.glob(diffPath + '/*.txt')
    
    for f in diffCurves:
        removeExt = f.split('.')
        fn = removeExt[0].split('\\')
        
        if fn[-1] is not 'Unfaulted_sim_exptDiff':
            faultedDiffQ, faultedDiff = pf.importFile(diffPath, fn[-1])
            fitDiff = np.subtract(unfaultedDiff, faultedDiff)
            
            
            saveName = fn[-1].split('_')
            with open(simPath + 'fitDiffCurves/' + saveName[0] + '_fitDiff.txt', 'w') as x:
                for (d) in fitDiff:
                    x.write('{0} \n'.format(d))
                x.close()
                
    return