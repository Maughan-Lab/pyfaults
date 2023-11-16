##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import os
import pandas as pd

#---------------------------------------------------------------------------------
# generates all supercell models within a given parameter space ------------------
#---------------------------------------------------------------------------------
def genSupercells(unitcell, nStacks, fltLayer, probList, sVecList, path):
    '''
    Parameters
    ----------
    unitcell
        Unitcell : base unit cell stack
    nStacks
        int : number of stacks per supercell
    fltLayer
        str : name of faulted layer
    probList
        list (float) : stacking fault probabilities to search
    sVecList
        list (nparray) : stacking vectors to search
    path
        str : CIF file save directory

    Returns
    -------
    df
        DataFrame : tabulated supercell data, includes the following as columns
            'Model' -- unique identifier
            'Stacking Vector' -- stacking vector in string format (TeX math mode)
            'Stacking Probability' -- stacking fault probability in str format (TeX math mode)
            'S_x' -- x-component of stacking vector
            'S_y' -- y-component of stacking vector
            'S_z' -- z-component of stacking vector
            'P' -- fault probability (0 to 1)
    '''
    import pyfaults as pf
    
    # create 'supercell_CIFs' folder in file directory
    if os.path.exists(path + 'supercell_CIFs/') == False:
        os.mkdir(path + 'supercell_CIFs/')
    
    df = pd.DataFrame()
    cellList = []
    probStrList, sVecStrList = pf.formatStr(probList, sVecList)
    
    # generate unfaulted supercell
    UF = pf.supercell.Supercell(unitcell, nStacks)
    cellList.append([UF, 'Unfaulted'])
    modelCol = ['Unfaulted']
    sVecTxtCol = ['None']
    probTxtCol = ['None']
    sxCol = [0]
    syCol = [0]
    szCol = [0]
    probCol = [0]
    
    # generate each faulted supercell
    for p in range(len(probList)):
        for s in range(len(sVecList)):
            FLT = pf.supercell.Supercell(unitcell, nStacks, fltLayer=fltLayer, 
                                         stackVec=sVecList[s], 
                                         stackProb=probList[p])
            cellTag = 'S' + str(s+1) + '_P' + str(int(probList[p]*100))
            cellList.append([FLT, cellTag])
            modelCol.append(cellTag)
            sVecTxtCol.append(sVecStrList[s])
            probTxtCol.append(probStrList[p])
            sxCol.append(sVecList[s][0])
            syCol.append(sVecList[s][1])
            szCol.append(sVecList[s][2])
            probCol.append(probList[p])
    
    # export CIF for each supercell
    for c in range(len(cellList)):
        pf.toCif((cellList[c][0]), path + 'supercell_CIFs/', cellList[c][1])
    
    # add entry for each supercell
    df['Model'] = modelCol
    df['Stacking Vector'] = sVecTxtCol
    df['Stacking Probability'] = probTxtCol
    df['S_x'] = sxCol
    df['S_y'] = syCol
    df['S_z'] = szCol
    df['P'] = probCol
    
    df.to_csv(path + 'supercell_CIFs/model_info.csv')
    
    return df               