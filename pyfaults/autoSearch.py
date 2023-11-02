########################################################################
# Author: Sinclair R. Combs
########################################################################

import numpy as np
import pandas as pd

#-----------------------------------------------------------------------
# automated parameter grid searching -----------------------------------
#-----------------------------------------------------------------------
def autoSearch(savePath, unitcell, nStacks, fltLayer, probList, sVecList,
               wl, maxTT, simPW, exptPath, exptName, peaksPath, peaksName):
    
    import pyfaults as pf
    
    # generate supercells ----------------------------------------------
    supercells = pf.genSupercells.genSupercells(unitcell, nStacks, fltLayer, 
                                                probList, sVecList, savePath)
    print('Finished generating supercell CIFs')
    
    # simulate XRD -----------------------------------------------------
    simData = pf.calcSims.calcSims(savePath + 'supercell_CIFs/', 'model_info',
                                   wl, maxTT, simPW, savePath)
    print('Finished simulating XRD patterns')
    
    # calculate difference curves --------------------------------------
    diffData = pf.calcDiffs.calcDiffs(exptPath, exptName, savePath + 'sims/',
                                      'sim_info', wl, maxTT)
    print('Finished calculating difference curves')
    
    # calculate fit differences ----------------------------------------
    fitDiffData = pf.calcFitDiffs.calcFitDiffs(savePath + 'sims/',
                                               savePath + 'diffCurves/', 
                                               'diffCurve_info')
    print('Finished calculating fit difference curves')
    
    # compare peak fits ------------------------------------------------
    compareFits = pf.compareFits.compareFits(savePath + 'fitDiffCurves/',
                                             'fitDiff_info', peaksPath, peaksName)
    print('Finished comparison of model fits')
    
    # display results --------------------------------------------------
    for i in compareFits.index:
        print(compareFits['Model'][i] + ' -- ' + compareFits['(hkl)'][i] +\
              '(' + compareFits['Peak Index (Q)'] + ' A^-1)' +\
              ' -- ' + compareFits['Faulted Model Fit'])
    
    return









        
        