##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import numpy as np
import os, glob
import pandas as pd

#---------------------------------------------------------------------------------
# generates simulated PXRD data from PyFaults input file -------------------------
#---------------------------------------------------------------------------------
def simulate(filePath):
    import pyfaults as pf

    structName = ''
    simType = ''
    gridSearch = ''
    latt = pf.lattice.Lattice(0,0,0,0,0,0)
    numLyrs = 0
    lyrNames = []
        
    f = open(filePath, 'r')
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
            
    #-----------------------------------------------------------------------------
    
    # generates unit cell and CIF
    unitcell = pf.unitcell.Unitcell(structName, lyrs, latt)
    unitcell.toCif('./')
    
    #-----------------------------------------------------------------------------

    fltLyr = ''
    prob = []
    sVec = []
    numStacks = 0
    supercells = []
    
    # creates folder to store supercell CIFs
    if os.path.exists('./supercells/') == False:
            os.mkdir('./supercells/')
            
    #-----------------------------------------------------------------------------
    gsProb = []
    gsSVec = []
    
    if gridSearch.startswith('Step'):
        pRange = []
        sxRange = []
        syRange = []
        
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
        
        gsProb, gsSVec = pf.gridSearch.stepGridSearch(pRange, sxRange, syRange)
        
    if gridSearch.startswith('Random'):
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
        
        gsProb, gsSVec = pf.gridSearch.randGridSearch(pRange, sxRange, syRange, numVec)
    
    #-----------------------------------------------------------------------------
    
    # displacement type SFs
    if simType.startswith('Displacement'):
        for i in range(len(inputFile)):
            
            # name of faulted layer
            if inputFile[i].startswith('FAULT LAYER:'):
                fltLyr = inputFile[i].split(':')[1]
                fltLyr = fltLyr.replace('\n', '')
            
            # number of stacks
            if inputFile[i].startswith('N:'):
                numStacks = int(inputFile[i].split(':')[1])
            
            # list of fault probabilities to simulate
            if gridSearch.startswith('None'):
                if inputFile[i].startswith('PROBABILITY:'):
                    removeVar = inputFile[i].split(':')[1]
                    pList = removeVar.split(',')
                    for j in pList:
                        prob.append(float(j))
            
            if gridSearch.startswith('Step'):
                prob = gsProb
            if gridSearch.startswith('Random'):
                prob = gsProb
            
            # list of displacement vectors to simulate
            if gridSearch.startswith('None'):
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
                            
            if gridSearch.startswith('Step'):
                sVec = gsSVec
            if gridSearch.startswith('Random'):
                sVec = gsSVec
        
        # generates faultless supercell
        noFlt = pf.supercell.Supercell(unitcell, numStacks)
        supercells.append(['Faultless', noFlt])
        
        # generates supercells for all combos of probabilities and vectors
        for p in range(len(prob)):
            for s in range(len(sVec)):
                newSupercell = pf.supercell.Supercell(unitcell, numStacks, fltLyr, sVec[s], prob[p])
                cellTag = 'S' + str(s+1) + '_P' + str(int(prob[p]*100))
                supercells.append([cellTag, newSupercell])
        
        # generates supercell CIFs
        for i in range(len(supercells)):
            pf.toCif(supercells[i][1], './supercells/', supercells[i][0])
            
    #-----------------------------------------------------------------------------
            
    # transition matrix type SFs
    if simType.startswith('Transition'):
        import random as r
        
        tmStartLyr = []
        tmNextLyr = []
        tmProb = []
        tmXvec = []
        tmYvec = []
        tmZvec = []
        
        for i in range(len(inputFile)):
        
            # number of stacks
            if inputFile[i].startswith('N:'):
                numStacks = int(inputFile[i].split(':')[1])
                
            # name of faulted layer
            if inputFile[i].startswith('FAULT LAYER'):
                fltLyr = inputFile[i].split(':')[1]
                fltLyr = fltLyr.replace('\n', '')
                
            # list of fault probabilities to simulate
            if inputFile[i].startswith('PROBABILITY'):
                removeVar = inputFile[i].split(':')[1]
                pList = removeVar.split(',')
                for j in pList:
                    prob.append(float(j))
                
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
            
            for stackProb in prob:
            
                seq = []
                seq.append(transMatrix['Start Layer'][0])
                
                numUCLyrs = len(lyrs)
                
                for n in range(numUCLyrs * numStacks):
                    p = (r.randint(0,100) * 0.01)
    
                    possLyrIndex = []
                    for i in range(len(transMatrix.index)):
                        lyrName = transMatrix['Start Layer'][i]
                        if lyrName == seq[-1]:
                            possLyrIndex.append(i)
                
                    for l in possLyrIndex:
                        isNext = False
                    
                        lyrProb = transMatrix['P'][l]
                        
                        if isinstance(lyrProb, str):
                            if lyrProb == 'P' or lyrProb == 'p':
                                if p <= stackProb:
                                    isNext = True
                                elif '-' in lyrProb:
                                    sub = lyrProb.split('-')
                                    if p > stackProb and p < float(sub[0]):
                                        isNext = True
                                elif float(lyrProb) == 1.0:
                                    isNext = True
                                
                        if isNext == True:
                            seq.append(transMatrix['Next Layer'][l])
                            break
                        
                newTMLyrs = []
                newTMLyrs.append(lyrs[0])
            
                for i in range(1, len(seq)):
                    lyrName = seq[i]
                    if lyrName == 'F':
                        lyrName = fltLyr
    
                    if i%2 == 0:
                        n = i/2
                    elif i%2 != 0:
                        n = (i+1)/2
    
                    adj = []
                    for j in range(len(transMatrix.index)):
                        if transMatrix['Start Layer'][j] == lyrName:
                            if transMatrix['Next Layer'][j] == seq[i+1]:
                                x = transMatrix['x'][j]
                                y = transMatrix['y'][j]
                                z = transMatrix['z'][j]
                                adj = [x, y, z*n]
                    
                    newLyrAtoms = []
                    for atom in lyrDict[lyrName]:
                        for a in atomDict:
                            if a==atom:
                                newAtom = pf.layerAtom.LayerAtom(lyrName, atom, atomDict[a][0],
                                                                 [float(atomDict[a][1]) + adj[0], 
                                                                  float(atomDict[a][2]) + adj[1], 
                                                                  float(atomDict[a][3]) + adj[2]], 
                                                                 float(atomDict[a][4]), 
                                                                 float(atomDict[a][5]), latt)
                                newLyrAtoms.append(newAtom)
                    newLyr = pf.layer.Layer(newLyrAtoms, latt, lyrName)
                    newTMLyrs.append(newLyr)
                    
                    cellTag = 'P' + str(int(prob[p]*100))
                    
                    cell = pf.unitcell.Unitcell(cellTag, newTMLyrs, latt)
                    supercell = pf.supercell.Supercell(cell, 1)
                    supercells.append([cellTag, supercell])
        
    #-----------------------------------------------------------------------------
    
    # intercalation type SFs
    if simType.startswith('Intercalation'):
        
        zAdj = 0
        for i in range(len(inputFile)):
            
            # adjustment to z-coordinate of faulted layer
            if inputFile[i].startswith('Z ADJ:'):
                zAdj = inputFile[i].split(':')[1]
                
            # name of faulted layer
            if inputFile[i].startswith('FAULT LAYER:'):
                fltLyr = inputFile[i].split(':')[1]
                fltLyr = fltLyr.replace('\n', '')
                
            # number of stacks
            if inputFile[i].startswith('N:'):
                numStacks = int(inputFile[i].split(':')[1])
            
            # intercalation layer
            if inputFile[i].startswith('INTLAYER:'):
                removeVar = inputFile[i].split(':')[1]
                intLayerAtoms = removeVar.split(',')
                intLayerAtoms[-1] = intLayerAtoms[-1].replace('\n', '')
                lyrDict['I'] = intLayerAtoms
                
            for a in lyrDict['I']:
                if inputFile[i].startswith(a):
                    removeVar = inputFile[i].split(':')[1]
                    atom = removeVar.split(',')
                    atom[-1] = atom[-1].replace('\n', '')
                    atomDict[a] = atom
                    
            lyrAtoms = []
            for atom in lyrDict['I']:
                for a in atomDict:
                    if a==atom:
                        newAtom = pf.layerAtom.LayerAtom(layer, atom, atomDict[a][0], 
                                                         [float(atomDict[a][1]), 
                                                          float(atomDict[a][2]), 
                                                          float(atomDict[a][3])], 
                                                         float(atomDict[a][4]), 
                                                         float(atomDict[a][5]), latt)
                        lyrAtoms.append(newAtom)
            newLayer = pf.layer.Layer(lyrAtoms, latt, 'I')
            lyrs.append(newLayer)
            
            # list of fault probabilities to simulate
            if inputFile[i].startswith('PROBABILITY:'):
                removeVar = inputFile[i].split(':')[1]
                pList = removeVar.split(',')
                for j in pList:
                    prob.append(float(j))
                
            # generates faultless supercell
            noFlt = pf.supercell.Supercell(unitcell, numStacks)
            supercells.append(['Faultless', noFlt])
            
            # generates supercells with intercalation layer added
            for p in range(len(prob)):
                newSupercell = pf.supercellAdjZ.SupercellAdjZ(unitcell, numStacks, lyrs[-1], fltLyr, zAdj, prob)
                cellTag = 'P' + str(int(prob[p]*100))
                supercells.append([cellTag, newSupercell])
                
            # generates supercell CIFs
            for i in range(len(supercells)):
                pf.toCif(supercells[i][1], './supercells/', supercells[i][0])
    
    #-----------------------------------------------------------------------------
                
    wl = 0.0
    maxTT = 0.0
    pw = 0.0
    
    # creates folder to store simulation data
    if os.path.exists('./simulations') == False:
        os.mkdir('./simulations')
    
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

    fileList = glob.glob('./supercells/*')
    for f in range(len(fileList)):
        fileList[f] = fileList[f].replace('.cif', '')
        fileList[f] = fileList[f].replace('./supercells\\', '')
    
    # PXRD simulation for each supercell
    for f in fileList:
        q, ints = pf.simXRD.fullSim('./supercells/', f, wl, maxTT, pw=pw, savePath='./simulations/')
        
    #-----------------------------------------------------------------------------
    
    # calculated R^2 values against experimental data if provided
    exptPath = ''
    exptFN = ''
    exptWL = 0
    
    for i in range(len(inputFile)):
        
        # file directory
        if inputFile[i].startswith('EXPT PATH'):
            exptPath = inputFile[i].split(':')[1]
            exptPath = exptPath.replace('\n', '')
        
        # experimental data file name
        if inputFile[i].startswith('EXPT FILE NAME'):
            exptFN = inputFile[i].split(':')[1]
            exptFN = exptFN.replace('\n', '')
        
        # instrument wavelength
        if inputFile[i].startswith('EXPT WAVELENGTH'):
            exptWL = float(inputFile[i].split(':')[1])
    
    # calculates and prints to text file
    if exptPath != '':
        pf.analyze.simR2vals('./simulations/', exptPath, exptFN, exptWL, maxTT)

    return