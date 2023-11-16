##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import os
import glob
import numpy as np

#---------------------------------------------------------------------------------
# calculates difference of fit residuals for all models in directory -------------
#---------------------------------------------------------------------------------
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
    
    # create 'fitDiffCurves' folder in file directory
    if os.path.exists(simPath + 'fitDiffCurves/') == False:
        os.mkdir(simPath + 'fitDiffCurves/')
    
    # import difference curves
    unfaultedDiffQ, unfaultedDiff = pf.importFile(diffPath, 'Unfaulted_sim_exptDiff')
    diffCurves = glob.glob(diffPath + '/*.txt')
    
    for f in diffCurves:
        removeExt = f.split('.')
        fn = removeExt[0].split('\\')
        
        # calculate difference of fit residuals
        if fn[-1] != 'Unfaulted_sim_exptDiff':
            faultedDiffQ, faultedDiff = pf.importFile(diffPath, fn[-1])
            fitDiff = np.subtract(unfaultedDiff, faultedDiff)
            
            # save calculated curve
            saveName = fn[-1].split('_exptDiff')
            with open(simPath + 'fitDiffCurves/' + saveName[0] + '_fitDiff.txt', 'w') as x:
                for (q, ints) in zip(faultedDiffQ, fitDiff):
                    x.write('{0} {1} \n'.format(q, ints))
                x.close()
                
    return