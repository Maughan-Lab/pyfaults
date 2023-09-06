#########################################################################################
# pyfaults.autoSearch
# Author: Sinclair R. Combs
#########################################################################################

import numpy as np
import pandas as pd
import re
from matplotlib.pyplot import rc
rc("text", usetex=True)
rc("font", **{"family":"sans-serif","sans-serif":["Helvetica"]},size="14")
rc("text.latex",preamble=r"\usepackage{sfmath}")


#----------------------------------------------------------------------------------------
def formatStr(probList, sVecList):
    
    probStrList = []
    for p in range(len(probList)):
        probStr = re.sub("x", str(int(probList[p]*100)), r"$P = x \%$")
        probStrList.append(probStr)
        
    sVecStrList = []
    for s in range(len(sVecList)):
        txt = r"$\vec{S}_n = \left[ x, y, z \right]$"
        sVecStr = re.sub("n", str(s+1), txt)
        strVars = ["x", "y", "z"]
        for v in range(len(strVars)):
            txtSub = re.sub(strVars[v], sVecList[s][v], sVecStr)
            sVecStr = txtSub
        sVecStrList.append(sVecStr)
        
    return probStrList, sVecStrList


#----------------------------------------------------------------------------------------
def genSupercells(unitcell, nStacks, fltLayer, probList, sVecList, path):
    from pyfaults.supercell import Supercell
    from pyfaults import toCif

    cellList = []
    
    probStrList, sVecStrList = formatStr(probList, sVecList)
    
    df = pd.DataFrame({"Model": [],"Stacking Vector": [], 
                       "Stacking Probability" : [],
                       "S_x": [], "S_y": [], "S_z": [], "P": []})
    
    UF = Supercell(unitcell, nStacks)
    cellList.append([UF, "Unfaulted"])
    
    df.loc[len(df)] = {"Model": "Unfaulted", 
                       "Stacking Vector": r"$\vec{S}_0 = \left[0, 0, 0 \right]$",
                       "Stacking Probability" : r"$P = 0 \%$",
                       "S_x": 0, "S_y": 0, "S_z": 0, "P": 0}
    
    for p in range(len(probList)):
        for s in range(len(sVecList)):
            FLT = Supercell(unitcell, nStacks, fltLayer=fltLayer, 
                            stackVec=sVecList[s], stackProb=probList[p])
            
            cellTag = "S" + str(s+1) + "_P" + str(int(probList[p]*100))
            cellList.append([FLT, cellTag])
            
            sx = sVecList[s][0]
            sy = sVecList[s][1]
            sz = sVecList[s][2]
            
            df.loc[len(df)] = {"Model": cellTag, 
                               "Stacking Vector": sVecStrList[s],
                               "Stacking Probability" : probStrList[p],
                               "S_x": sx, "S_y": sy, "S_z": sz, 
                               "P": probList[p]}
            
    for c in range(len(cellList)):
        toCif((cellList[c][0]), path, cellList[c][1])
            
    return df


#----------------------------------------------------------------------------------------
def calcSims(df, path, wl, maxTT, pw):
    import pyfaults.simXRD as xs
    
    simQ = []
    simInts = []
    for i in df.index:
        name = df["Model"][i]
        Q, ints = xs.fullSim(path, name, wl, maxTT, pw=pw)
        
        simQ.append(Q)
        normInts = xs.norm(ints)
        simInts.append(normInts)
    
    df.loc[:, ["Simulated Q"]] = simQ
    df.loc[:, ["Simulated Intensity"]] = simInts
    
    simDF = df
    
    return simDF


#----------------------------------------------------------------------------------------
def normalizeExpt(expt):
    import pyfaults.simXRD as xs
    
    exptNorm = [expt[0], xs.norm(expt[1])]
        
    return exptNorm


#----------------------------------------------------------------------------------------
def calcDiffs(path, simDF, expt):
    import pyfaults.simXRD as xs
    
    exptDiffs = []
    for i in simDF.index:
        name = simDF["Model"][i]
        qVals = simDF["Simulated Q"][i]
        intVals = simDF["Simulated Intensity"][i]
        
        diffQ, diffInts = xs.diffCurve(expt[0], qVals, expt[1], intVals)
        
        exptDiffs.append(diffInts)
        
        with open(path + "diffCurves/" + name + "_exptDiff.txt", "w") as f:
            for (diffInts) in zip(diffInts):
                f.write("{0}\n".format(diffInts))
        f.close() 
    
    simDF.loc[:, ["Expt vs. Model Difference"]] = exptDiffs
    exptDiffDF = simDF
    
    return exptDiffDF


#----------------------------------------------------------------------------------------
def calcFitDiffs(path, exptDiffDF):
    import pyfaults.simXRD as xs
    
    Q = exptDiffDF["Simulated Q"][0]
    UFdiff = exptDiffDF["Expt vs. Model Difference"][0]
    
    fitDiffs = []
    for i in range(1, exptDiffDF.index):
        name = exptDiffDF["Model"][i]
        modelDiff = exptDiffDF["Expt vs. Model Difference"][i]
        fitDiffQ, fitDiffInts = xs.diffCurve(Q, UFdiff, Q, modelDiff)
        
        fitDiffs.append(fitDiffInts)
        
        with open(path + "diffCurves/" + name + "_fitDiff.txt", "w") as f:
            for (fitDiffInts) in zip(fitDiffInts):
                f.write("{0}\n".format(fitDiffInts))
        f.close() 
    
    exptDiffDF.loc[:, ["Unfaulted vs. Faulted"]] = fitDiffs
    fitDiffDF = exptDiffDF
    
    return fitDiffDF


#----------------------------------------------------------------------------------------
def makePeakDF(labels, qVal, pw):
    df = pd.DataFrame({"Reflection": [],"Q": [], "Peak Range" : []})
    
    for i in range(len(labels)):
        qMin = qVal[i] - (pw/2)
        qMax = qVal[i] + (pw/2)
        qRange = [qMin, qMax]
        
        addRow = {"Reflection": labels[i], "Q": qVal[i], "Peak Range": qRange}
        df.loc[len(df)] = addRow
    
    return df


#----------------------------------------------------------------------------------------
def peakFitCompare(fitDiffDF, peakDF):
    
    Q = fitDiffDF["Simulated Q"][0]
    
    fitComp = []
    
    # for each model
    for i in fitDiffDF.index:
        name = fitDiffDF["Model"][i]
        
        compDF = pd.DataFrame({"Peak": [], "UF vs. FLT Fit Difference": [], 
                              "Sum": [], "Fit Result": []})
        
        fullFitDiff = fitDiffDF["Unfaulted vs. Faulted"][i]
        
        # for each peak
        for j in peakDF.index:
            peakLabel = peakDF["Reflection"][j]
            qRange = peakDF["Peak Range"][j]
            qMin = qRange[0]
            qMax = qRange[1]
            
            truncFitDiff = []
            for val in range(len(Q)):
                if Q[val] >= qMin and Q[val] <= qMax:
                    truncFitDiff.append(fullFitDiff[val])
                    
            diffSum = np.cumsum(truncFitDiff)
            if diffSum[-1] < 0:
                result = "Improved"
            elif diffSum[-1] > 0:
                result = "Worsened"
            elif diffSum[-1] == 0:
                result = "No Change"
                
            addRow = {"Peak": peakLabel, 
                      "UF vs. FLT Fit Difference": truncFitDiff,
                      "Sum": diffSum, "Fit Result": result}
            
            compDF.loc[len(compDF)] = addRow
            
        addEntry = [name, compDF]
        fitComp.append(addEntry)
        
    return fitComp


#----------------------------------------------------------------------------------------
def plotSearchResults(fitDiffList, peakQList, peakLabels, probList, sVecList, 
                      sVecLabels, xSpacing, wl):
    from pyfaults import plotXRD
    
    rows = len(fitDiffList)
    cols = len(peakQList)
    
    diffQ = []
    diffInts = []
    
    yMin = 0
    yMax = 0
    
    for i in range(rows):
        diffQ.append(fitDiffList[i][0])
        diffInts.append(fitDiffList[i][1])

        if np.max(fitDiffList[i][1]) > yMax:
            yMax = fitDiffList[i][1]
        if np.min(fitDiffList[i][0]) < yMin:
            yMin = fitDiffList[i][0]
            
    yLim = [yMin, yMax]
        
    xLims = []
    for i in range(cols):
        xLims.append([peakQList[i]-xSpacing, peakQList[i]+xSpacing])
    
    ax1 = plotXRD.fitCompare(rows, cols, diffQ, diffInts, xLims, yLim, wl, 
                             sVecLabels, peakLabels)
    
    return ax1


#----------------------------------------------------------------------------------------
def autoSearch(path, unitcell, expt, nStacks, fltLayer, probList, sVecList, 
               wl, maxTT, pw=0.0):
    
    cellList = genSupercells(unitcell, nStacks, fltLayer, probList, sVecList)
    
    simList = calcSims(path, wl, maxTT, pw)
    
    exptNorm, simListNorm = normalizeInts(expt, simList)
    
    exptDiffList = calcDiffs(path, simListNorm, exptNorm)
    
    fitDiffList = calcFitDiffs(path, exptDiffList)
    
    return fitDiffList

















        
        