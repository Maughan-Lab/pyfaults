#####################################################################################################################
# Author: Sinclair R. Combs
#####################################################################################################################

#--------------------------------------------------------------------------------------------------------------------
# automated parameter grid searching --------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
def autoSearch(path, unitcell, n, fltLayer, pList, sList, wl, maxTT, simPW, exptPath, exptName, brLabels, brQ, brPW):
    '''
    Parameters
    ----------
    path
        str : working directory 
    unitcell
        Unitcell : base unit cell used to generate supercells
    n
        int : number of stacks in supercells
    fltLayer
        str : name of layer to apply stacking vector to
    pList
        list (float) : list of fault probabilities to include in search
    sList
        list (nparray) : list of [x, y, z] stacking vectors to include in search
    wl
        float : instrument wavelength for PXRD simulations (A)
    maxTT
        float : maximum two theta value to calculate in PXRD simulations
    simPW
        float : artificial peak broadening term for PXRD simulations
    exptPath
        str : directory where experimental PXRD data file is stored
    exptName
        str : name of experimental PXRD data file, including file extension
    brLabels
        list (str) : list of (hkl) labels (strings) of broadened reflections
    brQ
        list (float) : list of Q (A^-1) indices of broadened reflections
    brPW
        float : additional Q (A^-1) range around indices for broadened reflections
        
    Returns
    -------
    supercells
        DataFrame : tabulated supercell model data, includes the following as columns
            'Model' -- unique identifier
            'Stacking Vector' -- stacking vector in string format (TeX math mode)
            'Stacking Probability' -- stacking fault probability in str format (TeX math mode)
            'S_x' -- x-component of stacking vector
            'S_y' -- y-component of stacking vector
            'S_z' -- z-component of stacking vector
            'P' -- fault probability (0 to 1)
    '''
    
    import pyfaults as pf
    
    # generate supercells ----------------------------------------------------------------
    supercells = pf.genSupercells.genSupercells(unitcell, n, fltLayer, pList, sList, path)
    print('Finished generating supercell CIFs')
    
    # simulate XRD -----------------------------------------------------------------------
    pf.calcSims.calcSims(path + 'supercell_CIFs/', 'model_info', wl, maxTT, simPW, path)
    print('Finished simulating XRD patterns')
    
    # calculate difference curves --------------------------------------------------------
    pf.calcDiffs.calcDiffs(path, exptName, path + 'sims/', wl, maxTT)
    print('Finished calculating difference curves')
    
    # calculate fit differences ----------------------------------------------------------
    pf.calcFitDiffs.calcFitDiffs(path + 'sims/', path + 'sims/diffCurves/')
    print('Finished calculating fit difference curves')
    
    # compare peak fits ------------------------------------------------------------------
    pf.compareFits.peakParams(brLabels, brQ, path, pw=brPW)
    pf.compareFits.compareFits(path + 'sims/fitDiffCurves/', path, 'peak_info')
    
    return supercells









        
        