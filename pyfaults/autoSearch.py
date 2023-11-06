########################################################################
# Author: Sinclair R. Combs
########################################################################


#-----------------------------------------------------------------------
# automated parameter grid searching -----------------------------------
#-----------------------------------------------------------------------
def autoSearch(savePath, unitcell, nStacks, fltLayer, probList, sVecList,
               wl, maxTT, simPW, exptPath, exptName, peaksCSVPath, peaksCSVName):
    
    import pyfaults as pf
    
    # generate supercells ----------------------------------------------
    supercells = pf.genSupercells.genSupercells(unitcell, nStacks, fltLayer, 
                                                probList, sVecList, savePath)
    print('Finished generating supercell CIFs')
    
    # simulate XRD -----------------------------------------------------
    pf.calcSims.calcSims(savePath + 'supercell_CIFs/', 'model_info', wl, maxTT, 
                         simPW, savePath)
    print('Finished simulating XRD patterns')
    
    # calculate difference curves --------------------------------------
    pf.calcDiffs.calcDiffs(exptPath, exptName, savePath + 'sims/', wl, maxTT)
    print('Finished calculating difference curves')
    
    # calculate fit differences ----------------------------------------
    pf.calcFitDiffs.calcFitDiffs(savePath + 'sims/', savePath + 'sims/diffCurves/')
    print('Finished calculating fit difference curves')
    
    # compare peak fits ------------------------------------------------
    pf.compareFits.compareFits(savePath + 'sims/fitDiffCurves/',
                               peaksCSVPath, peaksCSVName)
    
    return supercells









        
        