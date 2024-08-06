##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import numpy as np

def stepGridSearch(pRange, sxRange, syRange):
    '''
    Parameters
    ----------
    pRange (array_like) : list of minimum fault probability, maximum fault probability, and step size
    sxRange (array_like) : list of minimum stacking vector x-value, maximum stacking vector x-value, and step size
    syRange (array_like) : list of minimum stacking vector y-value, maximum stacking vector y-value, and step size
        
    Returns
    -------
    pList (array_like) : list of probabilities
    sList (array_like) : list of vectors
    '''
    
    # generate fault probabilities
    p = pRange[0]
    pList = []
    while p <= pRange[1]:
        pList.append(round(p, 3))
        p = p + pRange[2]
    
    # generate stacking vectors
    sx = sxRange[0]
    sy = syRange[0]
    sxList = []
    syList = []
    sList = []
    while sx <= sxRange[1]:
        sxList.append(round(sx, 5))
        sx = sx + sxRange[2]
    while sy <= syRange[1]:
        syList.append(round(sy, 5))
        sy = sy + syRange[2]
    for i in range(len(sxList)):
        for j in range(len(syList)):
            s = [sxList[i], syList[j], 0]
            sList.append(s)
            
    return np.array(pList), np.array(sList)

def randGridSearch(pRange, sxRange, syRange, numVec):
    '''
    Parameters
    -------
    pRange (array_like) : list of minimum fault probability, maximum fault probability, and step size
    sxRange (array_like) : list of minimum stacking vector x-value, maximum stacking vector x-value, and step size
    syRange (array_like) : list of minimum stacking vector y-value, maximum stacking vector y-value, and step size
    numVec (int) : number of randomized stacking vectors to generate
        
    Returns
    -------
    pList (array_like) : list of probabilities
    sList (array_like) : list of vectors
    '''
    import random
    
    # generate fault probabilities
    p = pRange[0]
    pList = []
    while p <= pRange[1]:
        pList.append(round(p, 3))
        p = p + pRange[2]
    
    # generate stacking vectors
    sList = []
    for i in range(numVec):
        sx = random.randrange(sxRange[0], sxRange[1])
        sy = random.randrange(syRange[0], syRange[1])
        sList.append(sx, sy, 0)
    
    return np.array(pList), np.array(sList)