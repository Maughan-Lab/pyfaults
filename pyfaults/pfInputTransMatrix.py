import random as r
import pandas as pd

def pfInputTransMatrix(path, prob, lyrs, numStacks, fltLyr, lyrDict, atomDict, latt):
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
        for n in range(total):
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
        newTMLyrs.append(lyrs[0])
    
        for i in range(1, len(seq)-1):
            lyrName = seq[i]
            if lyrName == 'F':
                lyrName = fltLyr

            if i%2 == 0:
                n = i/2
            elif i%2 != 0:
                n = (i+1)/2

            test = transMatrix[(transMatrix['Start Layer'] == lyrName) & (transMatrix['Next Layer'] == seq[i+1])]
            adj = [test.iloc[0]['x'], test.iloc[0]['y'], test.iloc[0]['z']*n]
                        
            newLyrAtoms = []
            for a in lyrDict[lyrName]:
                newAtom = pf.layerAtom.LayerAtom(lyrName, a, atomDict[a][0],
                                                             [float(atomDict[a][1]) + adj[0], 
                                                              float(atomDict[a][2]) + adj[1], 
                                                              float(atomDict[a][3]) + adj[2]], 
                                                             float(atomDict[a][4]), 
                                                             float(atomDict[a][5]), latt)
                newLyrAtoms.append(newAtom) 
            newLyr = pf.layer.Layer(newLyrAtoms, latt, lyrName)
            newTMLyrs.append(newLyr)

        newTMLyrs.pop(0)
        cell = pf.unitcell.Unitcell(currP, newTMLyrs, latt)
        cells.append([cell, currP])
            
    return cells
