##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import os

# generates all supercell models within a given parameter space ----------
def genSupercells(unitcell, nStacks, fltLayer, probList, sVecList):
    '''
    Parameters
    ----------
    unitcell (Unitcell) : base unit cell
    nStacks (int) : number of unit cell stacks in supercell
    fltLayer (str) : faulted layer name
    probList (array_like) : list of stacking fault probabilities
    sVecList (array_like) : list of displacement vectors [x,y,z] in fractional coordinates
    '''
    import pyfaults as pf
    
    # create 'supercells' folder in working directory
    if os.path.exists('./supercells/') == False:
        os.mkdir('./supercells/')
    
    cellList = []
    
    # generate unfaulted supercell
    UF = pf.supercell.Supercell(unitcell, nStacks)
    cellList.append([UF, 'Unfaulted'])
    
    # generate each faulted supercell
    for p in range(len(probList)):
        for s in range(len(sVecList)):
            FLT = pf.supercell.Supercell(unitcell, nStacks, fltLayer=fltLayer, stackVec=sVecList[s], stackProb=probList[p])
            cellTag = 'S' + str(s+1) + '_P' + str(int(probList[p]*100))
            cellList.append([FLT, cellTag])
    
    # export CIF for each supercell
    for c in range(len(cellList)):
        pf.toCif((cellList[c][0]), './supercells/', cellList[c][1])
    
    return               
