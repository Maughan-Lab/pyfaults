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
        for i in range(1, len(seq)-1):
            lyrName = seq[i]

            counter = counter + 1
            if counter == numUCLyrs:
                nCount = nCount + 1
                counter = 0

            test = transMatrix[(transMatrix['Start Layer'] == lyrName) & (transMatrix['Next Layer'] == seq[i+1])]
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
