##################################################################################
# Author: Sinclair R. Combs
##################################################################################

# read gridSearch input and generate parameter space ----------
def pfInputGridSearch(gsData):
    import pyfaults as pf
    
    gsProb = []
    gsSVec = []
    
    gridSearch = gsData.loc[0, 'gridSearch']
    pRange = gsData.loc[0, 'pRange']
    sxRange = gsData.loc[0, 'sxRange']
    syRange = gsData.loc[0, 'syRange']
    numVec = gsData.loc[0, 'numVec']
    
    if gridSearch.iloc[0] == 'Step':
        gsProb, gsSVec = pf.gridSearch.stepGridSearch(pRange.iloc[0], sxRange.iloc[0], syRange.iloc[0])
        
    if gridSearch.iloc[0] == 'Random':
        gsProb, gsSVec = pf.gridSearch.randGridSearch(pRange.iloc[0], sxRange.iloc[0], syRange.iloc[0], numVec.iloc[0])
        
    return gsProb, gsSVec
