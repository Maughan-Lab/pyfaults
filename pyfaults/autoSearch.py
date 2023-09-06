#########################################################################################
# pyfaults.autoSearch
# Author: Sinclair R. Combs
#########################################################################################

def genSupercells(unitcell, nStacks, fltLayer, probList, sVecList):
    from pyfaults.supercell import Supercell
    
    cellList = []
    
    UF = Supercell(unitcell, nStacks)
    cellList.append([UF, "Unfaulted", 0, "UF"])
    
    for p in range(len(probList)):
        probTag = str(int(probList[p] * 100))
        for s in range(len(sVecList)):
            FLT = Supercell(unitcell, nStacks, fltLayer=fltLayer, 
                            stackVec=sVecList[s], stackProb=probList[p])
            
            sx = str(sVecList[s][0])
            sy = str(sVecList[s][1])
            sz = str(sVecList[s][2])
            
            sVecTag = "[" + sx + ", " + sy + ", " + sz + "]"
            cellTag = "S" + str(s+1) + "_P" + probTag
            
            cellList.append([FLT, sVecTag, probTag, cellTag])
            
    return cellList


def cellListToCIF(cellList, path):
    from pyfaults import toCif
    
    for cell in range(len(cellList)):
        toCif(cellList[cell][0], path, cellList[cell][3])


def calcSims(path, wl, maxTT, pw):
    import glob
    import pyfaults.simXRD as xs
    
    CIFfileList = glob.glob(path + "*.cif")
        
    simList = []
    for i in range(len(CIFfileList)):
        removeExt = CIFfileList[i].split(".")
        cellName = removeExt[0].split("\\")
        Q, ints = xs.fullSim(path, cellName[1], wl, maxTT, pw=pw)
        
        simList.append([Q, ints, cellName[1]])
        xs.saveSim(path, cellName[1] + "_sim", Q, ints)
    
    return simList


def calcDiffs(path, simList, expt):
    import pyfaults.simXRD as xs
    
    exptDiffList = []
    for i in range(len(simList)):
        diffQ, diffInts = xs.diffCurve(expt[0], simList[i][0], expt[1], simList[i][1])
        
        exptDiffList.append([diffQ, diffInts, simList[i][2]])
        xs.saveDiff(path, "expt_" + simList[i][2] + "_diff", 
                    diffQ, diffInts)
    
    return exptDiffList


def calcFitDiffs(path, exptDiffList):
    import pyfaults.simXRD as xs
    
    UFdiffQ = exptDiffList[0][0]
    UFdiffInts = exptDiffList[0][1]
    
    fitDiffList = []
    for i in range(1, len(exptDiffList)):
        fitDiffQ, fitDiffInts = xs.diffCurve(UFdiffQ, exptDiffList[i][0], 
                                             UFdiffInts, exptDiffList[i][1])
        
        fitDiffList.append([fitDiffQ, fitDiffInts, exptDiffList[i][2]])
        xs.saveDiff(path, exptDiffList[i][2] + "_fitDiff", 
                    fitDiffQ, fitDiffInts)
    
    return fitDiffList


def autoSearch(path, unitcell, expt, nStacks, fltLayer, probList, sVecList, 
               wl, maxTT, pw=0.0):
    
    cellList = genSupercells(unitcell, nStacks, fltLayer, probList, sVecList)
    
    cellListToCIF(cellList, path)
    
    simList = calcSims(path, wl, maxTT, pw)
    
    exptDiffList = calcDiffs(path, simList, expt)
    
    fitDiffList = calcFitDiffs(path, exptDiffList)
    
    return fitDiffList

















        
        