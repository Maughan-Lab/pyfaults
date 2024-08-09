##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import pandas as pd

# read PyFaults input file ----------
def pfInput(path):
    '''
    Parameters
    ----------
    path (str) : file path where input file is stored
    '''
    
    import pyfaults as pf
    
    # initialize unitcell variables
    structName = ''
    simType = ''
    gridSearch = 'None'
    latt = pf.lattice.Lattice(0,0,0,0,0,0)
    numLyrs = 0
    lyrNames = []
    
    # read input file
    f = open(path, 'r')
    inputFile = f.readlines()
    
    for i in range(len(inputFile)):
        # structure name
        if inputFile[i].startswith('NAME'):
            structName = inputFile[i].split(':')[1]
            structName = structName.replace('\n', '')
        
        # type of stacking fault system: Displacement, Transition Matrix, or Intercalation
        if inputFile[i].startswith('TYPE'):
            simType = inputFile[i].split(':')[1]
            simType = simType.replace('\n', '')
            
        # grid search function: None, Step, or Random
        if inputFile[i].startswith('GRID SEARCH'):
            gridSearch = inputFile[i].split(':')[1]
            gridSearch = gridSearch.replace('\n', '')
        
        # lattice parameters
        if inputFile[i].startswith('LATTICE'):
            removeVar = inputFile[i].split(':')[1]
            strLatt= removeVar.split(',')
            strLatt[5] = strLatt[5].replace('\n', '')
            latt = pf.lattice.Lattice(a=float(strLatt[0]), 
                                      b=float(strLatt[1]), 
                                      c=float(strLatt[2]), 
                                      alpha=float(strLatt[3]), 
                                      beta=float(strLatt[4]), 
                                      gamma=float(strLatt[5]))
        
        # number of unique layers defined
        if inputFile[i].startswith('NUM LAYERS'):
            numLyrs = inputFile[i].split(':')[1]
            numLyrs = numLyrs.replace('\n', '')
      
    # generates layer name list (L1, L2, L3, etc.)
    for i in range(int(numLyrs)):
        lyrNames.append('L' + str(i+1))
    
    # generates dictionary of atoms within a layer
    lyrDict = {}
    for l in lyrNames:
        for i in range(len(inputFile)):
            if inputFile[i].startswith(l):
                removeVar = inputFile[i].split(':')[1]
                atomNames = removeVar.split(',')
                atomNames[-1] = atomNames[-1].replace('\n', '')
                lyrDict[l] = atomNames
    
    # generates dictionary of atomic parameters
    atomDict = {}
    for l in lyrNames:
        for a in lyrDict[l]:
            for i in range(len(inputFile)):
                if inputFile[i].startswith(a):
                    removeVar = inputFile[i].split(':')[1]
                    atom = removeVar.split(',')
                    atom[-1] = atom[-1].replace('\n', '')
                    atomDict[a] = atom
    
    # generates LayerAtom and Layer objects
    lyrs = []
    for layer in lyrNames:
        lyrAtoms = []
        for atom in lyrDict[layer]:
            for a in atomDict:
                if a==atom:
                    newAtom = pf.layerAtom.LayerAtom(layer, atom, atomDict[a][0], 
                                                     [float(atomDict[a][1]), 
                                                      float(atomDict[a][2]), 
                                                      float(atomDict[a][3])], 
                                                     float(atomDict[a][4]), 
                                                     float(atomDict[a][5]), latt)
                    lyrAtoms.append(newAtom)
        newLayer = pf.layer.Layer(lyrAtoms, latt, layer)
        lyrs.append(newLayer)
    
    # generates chemically identical layers in different positions in the unit cell
    copyLayers = []
    for i in range(len(inputFile)):
        if inputFile[i].startswith('COPY LAYER:'):
            removeVar = inputFile[i].split(':')[1]
            copyLyr = removeVar.split(',')
            lyrNames.append(copyLyr[0])
            adj = [float(copyLyr[2]), float(copyLyr[3]), float(copyLyr[4])]

            lyrDict[copyLyr[0]] = lyrDict[copyLyr[1]]

            copyLyrAtoms = []
            for atom in lyrDict[copyLyr[0]]:
                for a in atomDict:
                    if a==atom:
                        newAtom = pf.layerAtom.LayerAtom(copyLyr[0], atom, atomDict[a][0], 
                                                     [float(atomDict[a][1]) + adj[0], 
                                                      float(atomDict[a][2]) + adj[1], 
                                                      float(atomDict[a][3]) + adj[2]], 
                                                     float(atomDict[a][4]), 
                                                     float(atomDict[a][5]), latt)
                        copyLyrAtoms.append(newAtom)
            newCopyLayer = pf.layer.Layer(copyLyrAtoms, latt, copyLyr[0])
            lyrs.append(newCopyLayer)
            copyLayers.append(newCopyLayer)
            
    # generates unit cell
    unitcell = pf.unitcell.Unitcell(structName, lyrs, latt)
    
    # initialize gridSearch variables
    pRange = []
    sxRange = []
    syRange = []
    numVec = 0
    
    for i in range(len(inputFile)):
        if inputFile[i].startswith('P RANGE'):
            removeVar = inputFile[i].split(':')[1]
            pList = removeVar.split(',')
            pList[-1] = pList[-1].replace('\n', '')
            for j in pList:
                pRange.append(float(j))
            
            
        if inputFile[i].startswith('SX RANGE'):
            removeVar = inputFile[i].split(':')[1]
            sxList = removeVar.split(',')
            sxList[-1] = sxList[-1].replace('\n', '')
            for j in sxList:
                sxRange.append(float(j))
            
        if inputFile[i].startswith('SY RANGE'):
            removeVar = inputFile[i].split(':')[1]
            syList = removeVar.split(',')
            syList[-1] = syList[-1].replace('\n', '')
            for j in syList:
                syRange.append(float(j))
                
        if inputFile[i].startswith('NUM VECTORS'):
            numVec = inputFile[i].split(':')[1]
            numVec = numVec.replace('\n', '')
    
    # initialize supercell variables
    fltLyr = ''
    prob = []
    sVec = []
    numStacks = 0
    zAdj = 0
    
    for i in range(len(inputFile)):
        
        if inputFile[i].startswith('FAULT LAYER:'):
            fltLyr = inputFile[i].split(':')[1]
            fltLyr = fltLyr.replace('\n', '')
            
        if inputFile[i].startswith('PROBABILITY:'):
            removeVar = inputFile[i].split(':')[1]
            pList = removeVar.split(',')
            for j in pList:
                prob.append(float(j))
                
        if inputFile[i].startswith('VECTOR:'):
            removeVar = inputFile[i].split(':')[1]
            sList = removeVar.split(';')
            sList[-1] = sList[-1].replace('\n', '')
            for vec in sList:
                vec = vec[1:-1].split(',')
                arrayVec = []
                for j in vec:
                    j = float(j)
                    arrayVec.append(j)
                    sVec.append(arrayVec)
                    
        if inputFile[i].startswith('N:'):
            numStacks = int(inputFile[i].split(':')[1])
            
        if inputFile[i].startswith('Z ADJ:'):
            zAdj = inputFile[i].split(':')[1]
            
    # initialize simulation variables
    wl = 0.0
    maxTT = 0.0
    pw = 0.0
    
    for i in range(len(inputFile)):
        
        # simulation wavelength
        if inputFile[i].startswith('WAVELENGTH'):
            wl = float(inputFile[i].split(':')[1])
        
        # simulation maximum two theta in degrees
        if inputFile[i].startswith('MAX TWO THETA'):
            maxTT = float(inputFile[i].split(':')[1])
            
        # artificial broadening term
        if inputFile[i].startswith('BROADENING'):
            pw = float(inputFile[i].split(':')[1])
    
    
    ucData = [structName, simType, latt, numLyrs, lyrNames, lyrDict, atomDict, lyrs]
    ucCols = ['structName', 'simType', 'latt', 'numLyrs', 'lyrNames', 'lyrDict', 'atomDict', 'lyrs']
    ucDF = pd.DataFrame(data=[ucData], columns=[ucCols])
    
    gsData = [gridSearch, pRange, sxRange, syRange, numVec]
    gsCols = ['gridSearch', 'pRange', 'sxRange', 'syRange', 'numVec']  
    gsDF = pd.DataFrame(data=[gsData], columns=[gsCols])
    
    scData = [fltLyr, prob, sVec, numStacks, zAdj]
    scCols = ['fltLyr', 'prob', 'sVec', 'numStacks', 'zAdj']
    scDF = pd.DataFrame(data=[scData], columns=[scCols])
    
    simData = [wl, maxTT, pw]
    simCols = ['wl', 'maxTT', 'pw']
    simDF = pd.DataFrame(data=[simData], columns=[simCols])
            
    return unitcell, ucDF, gsDF, scDF, simDF
