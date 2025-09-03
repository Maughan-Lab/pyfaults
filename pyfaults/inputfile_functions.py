""" 
inputfile_functions.py

Module containing functions for utilizing a Pyfaults-style input file

pfInput --> 
pfInputGridSearch -->
pfInputTransMatrix --> 
"""

#---------- import packages ----------
import pandas as pd
import os



#-------------------------------------
#---------- FUNCTION: toCif ----------
#-------------------------------------
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
        # get structure name
        if inputFile[i].startswith('NAME'):
            structName = inputFile[i].split(':')[1]
            structName = structName.replace('\n', '')
        
        # get type of stacking fault system
        # must be Displacement, Transition Matrix, or Intercalation
        if inputFile[i].startswith('TYPE'):
            simType = inputFile[i].split(':')[1]
            simType = simType.replace('\n', '')
            
        # get grid search function
        # must be None, Step, or Random
        if inputFile[i].startswith('GRID SEARCH'):
            gridSearch = inputFile[i].split(':')[1]
            gridSearch = gridSearch.replace('\n', '')
        
        # get lattice parameters
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
        
        # get number of unique layers
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
    
    # generate supercells
    if simType == 'Displacement':
        pf.genSupercells.genSupercells(unitcell, numStacks, fltLyr, prob, sVec)
        
    elif simType == 'Transition Matrix':
        tmCells = pf.pfInputTransMatrix.pfInputTransMatrix(path, prob, lyrs, numStacks, fltLyr, lyrDict, atomDict, latt)
        
        if os.path.exists('./supercells/') == False:
            os.mkdir('./supercells/')
        
        for i in range(len(tmCells)):
            pf.toCif(tmCells[i][0], './supercells/', tmCells[i][1])
            
    return unitcell, ucDF, gsDF, scDF, simDF

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

def pfInputTransMatrix(path, prob, lyrs, numStacks, fltLyr, lyrDict, atomDict, latt):
    
    import pandas as pd
    import numpy as np
    import random as r
    import pyfaults as pf
    
    transMatrix = []
    
    f = open(path, 'r')
    inputFile = f.readlines()
    
    tmStartLyr = []
    tmNextLyr = []
    tmProb = []
    tmXvec = []  
    tmYvec = []
    tmZvec = []
    
    cells = []
    
    for i in range(len(inputFile)):
        if inputFile[i].startswith('TM:'):
            removeVar = inputFile[i].split(':')[1]
            transParams = removeVar.split(',')
            tmLyrs = transParams[0].split('-')
            
            tmStartLyr.append(tmLyrs[0])
            tmNextLyr.append(tmLyrs[1])
            tmProb.append(transParams[1])
            tmXvec.append(float(transParams[2]))
            tmYvec.append(float(transParams[3]))
            tmZvec.append(float(transParams[4]))
            
    d = {'Start Layer': tmStartLyr, 'Next Layer': tmNextLyr, 'P': tmProb,
          'x': tmXvec, 'y': tmYvec, 'z': tmZvec}
    transMatrix = pd.DataFrame(data=d)

    newLatt = pf.lattice.Lattice(latt.a, latt.b, (latt.c * numStacks), latt.alpha, latt.beta, latt.gamma)
    
    for stackProb in prob:
        seq = []
        seq.append(transMatrix.iloc[0]['Start Layer'])
        
        numUCLyrs = len(lyrs)
        
        # replace probability variables
        numProb = []
        numProbMin = []
        numProbMax = []
        for i in range(len(transMatrix.index)):
            lyrProb = transMatrix['P'][i]
            if isinstance(lyrProb, str):
                newProb = lyrProb
                newProb = newProb.replace('P', str(stackProb))
                newProb = newProb.replace('p', str(stackProb))

                if '-' in newProb:
                    sub = newProb.split('-')
                    newProb = float(sub[0])-float(sub[1])
                    numProb.append(float(newProb))
                    numProbMin.append(float(newProb))
                    numProbMax.append(float(sub[0]))
                else:
                    numProb.append(float(newProb))
                    numProbMin.append(0.0)
                    numProbMax.append(float(newProb))   
            else:
                numProb.append(float(newProb))
                numProbMin.append(0.0)
                numProbMax.append(newProb)
                
        currP = 'P' + str(int(stackProb*100))
        transMatrix[currP] = numProb
        transMatrix[currP + '_min'] = numProbMin
        transMatrix[currP + '_max'] = numProbMax

        total = int(numUCLyrs * numStacks)
        while len(seq) < total:
            p = (r.randint(0,100) * 0.01)

            possibleTrans = []
            for i in range(len(transMatrix.index)):
                if transMatrix['Start Layer'][i] == seq[-1]:
                    possibleTrans.append(i)

            nextLyr = []
            for i in possibleTrans:
                if transMatrix[currP][i] != 0.0:
                    min = transMatrix[currP + '_min'][i]
                    max = transMatrix[currP + '_max'][i]
                
                    if transMatrix[currP][i] == 1.0:
                        nextLyr.append(transMatrix['Next Layer'][i])
                    elif p >= min and p < max:
                        nextLyr.append(transMatrix['Next Layer'][i])
                    else:
                        nextLyr.append('None')
                else:
                    nextLyr.append('None')

            for i in nextLyr:
                if i != 'None':
                    seq.append(i)
                    
        newTMLyrs = []
        startLyrAtoms = []
        for a in lyrDict[lyrs[0].layerName]:
            newAtom = pf.layerAtom.LayerAtom(lyrs[0].layerName, a, atomDict[a][0],
                                                             [float(atomDict[a][1]), 
                                                              float(atomDict[a][2]), 
                                                              (float(atomDict[a][3])/numStacks)], 
                                                             float(atomDict[a][4]), 
                                                             float(atomDict[a][5]), newLatt)
            startLyrAtoms.append(newAtom)
        newStartLyr = pf.layer.Layer(startLyrAtoms, newLatt, lyrs[0].layerName)
        newTMLyrs.append(newStartLyr)

        counter = 0
        nCount = 0
        for i in range(1, len(seq)):
            lyrName = seq[i]

            counter = counter + 1
            if counter == numUCLyrs:
                nCount = nCount + 1
                counter = 0
                
            if i < len(seq)-1:
                test = transMatrix[(transMatrix['Start Layer'] == lyrName) & (transMatrix['Next Layer'] == seq[i+1])]
            elif i == len(seq):
                test = transMatrix[(transMatrix['Start Layer'] == lyrName) & (transMatrix['Next Layer'] == seq[0])]
            adj = [test.iloc[0]['x'], test.iloc[0]['y'], test.iloc[0]['z']]
            
            newLyrAtoms = []
            if lyrName == 'F':
                for a in lyrDict[fltLyr]:
                    newX = float(atomDict[a][1]) + adj[0]
                    newY = float(atomDict[a][2]) + adj[1]
                    newZ = np.round(float((float(atomDict[a][3])+nCount)/numStacks + adj[2]), 2)
                    
                    newAtom = pf.layerAtom.LayerAtom(lyrName, a, atomDict[a][0], [newX, newY, newZ], float(atomDict[a][4]), 
                                                     float(atomDict[a][5]), newLatt)
                    newLyrAtoms.append(newAtom)  
                    
                newLyr = pf.layer.Layer(newLyrAtoms, newLatt, lyrName)
                newTMLyrs.append(newLyr)
                
            else:    
                for a in lyrDict[lyrName]:
                    newX = float(atomDict[a][1]) + adj[0]
                    newY = float(atomDict[a][2]) + adj[1]
                    newZ = np.round(float((float(atomDict[a][3])+nCount)/numStacks + adj[2]), 2)
                        
                    newAtom = pf.layerAtom.LayerAtom(lyrName, a, atomDict[a][0], [newX, newY, newZ], float(atomDict[a][4]), 
                                                     float(atomDict[a][5]), newLatt)
                    newLyrAtoms.append(newAtom)     
                
                newLyr = pf.layer.Layer(newLyrAtoms, newLatt, lyrName)
                newTMLyrs.append(newLyr)

        cell = pf.unitcell.Unitcell(currP, newTMLyrs, newLatt)
        cells.append([cell, currP])

    return cells