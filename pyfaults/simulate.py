##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import numpy as np
import os, glob

#---------------------------------------------------------------------------------
# generates simulated PXRD data from PyFaults input file -------------------------
#---------------------------------------------------------------------------------
def simulate(filePath):
    import pyfaults as pf

    structName = ''
    simType = ''
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
        
        # type of stacking fault system: displacement, transition matrix, or intercalation
        if inputFile[i].startswith('TYPE'):
            simType = inputFile[i].split(':')[1]
        
        # lattice parameters
        if inputFile[i].startswith('LATTICE'):
            removeVar = inputFile[i].split(':')[1]
            strLatt= removeVar.split(',')
            latt = pf.lattice.Lattice(a=float(strLatt[0]), 
                                      b=float(strLatt[1]), 
                                      c=float(strLatt[2]), 
                                      alpha=float(strLatt[3]), 
                                      beta=float(strLatt[4]), 
                                      gamma=float(strLatt[5]))
        
        # number of unique layers defined
        if inputFile[i].startswith('NUM LAYERS'):
            numLyrs = inputFile[i].split(':')[1]
      
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
                lyrDict[l] = atomNames
    
    # generates dictionary of atomic parameters
    atomDict = {}
    for l in lyrNames:
        for a in lyrDict[l]:
            for i in range(len(inputFile)):
                if inputFile[i].startswith(a):
                    removeVar = inputFile[i].split(':')[1]
                    atom = removeVar.split(',')
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
                        newAtom = pf.layerAtom.LayerAtom(layer, atom, atomDict[a][0], 
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
    
    # displacement type SFs
    if simType.startswith('Displacement'):
        for i in range(len(inputFile)):
            
            # name of faulted layer
            if inputFile[i].startswith('FAULT LAYER:'):
                fltLyr = inputFile[i].split(':')[1]
            
            # number of stacks
            if inputFile[i].startswith('N:'):
                numStacks = int(inputFile[i].split(':')[1])
            
            # list of fault probabilities to simulate
            if inputFile[i].startswith('PROBABILITY:'):
                removeVar = inputFile[i].split(':')[1]
                pList = removeVar.split(',')
                for j in pList:
                    prob.append(float(j))
            
            # list of displacement vectors to simulate
            if inputFile[i].startswith('VECTOR:'):
                removeVar = inputFile[i].split(':')[1]
                sList = removeVar.split(';')
                for j in sList:
                    j.replace('[', '')
                    j.replace(']', '')
                    sVec.append(np.array(j))
        
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
    # if simType.startswith('Transition'):
        
    #-----------------------------------------------------------------------------
    
    # intercalation type SFs
    # if simType.startswith('Intercalation'):
    
    #-----------------------------------------------------------------------------
                
    wl = 0.0
    maxTT = 0.0
    pw = 0.0
    
    # creates folder to store simulation data
    if os.path.exists('./simulations') == False:
        os.mkdir('./simulations')
    
    for i in range(len(inputFile)):
        
        # simulation wavelength
        if inputFile[i].startswith('WAVELENGTH:'):
            wl = float(inputFile[i].split(':')[1])
        
        # simulation maximum two theta in degrees
        if inputFile[i].startswith('MAX TWO THETA:'):
            maxTT = float(inputFile[i].split(':')[1])
            
        # artificial broadening term
        if inputFile[i].startswith('BROADENING:'):
            pw = float(inputFile[i].split(':')[1])

    fileList = glob.glob('./supercells/*')
    for f in range(len(fileList)):
        fileList[f] = fileList[f].replace('.cif', '')
        fileList[f] = fileList[f].replace('./supercells\\', '')
    
    # PXRD simulation for each supercell
    for f in fileList:
        q, ints = pf.simXRD.fullSim('./supercells/', f, wl, maxTT, pw=pw, savePath='./simulations/')

    return