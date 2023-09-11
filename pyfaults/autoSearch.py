#########################################################################################
# pyfaults.autoSearch
# Author: Sinclair R. Combs
#########################################################################################

import numpy as np
import copy as cp
import pandas as pd
import re
import os

import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
rc("text", usetex=True)
rc("font", **{"family":"sans-serif","sans-serif":["Helvetica"]},size="14")
rc("text.latex",preamble=r"\usepackage{sfmath}")

#----------------------------------------------------------------------------------------
def importExpt(path, fn, wl, maxTT):
    from pyfaults import importSim, tt_to_q
    import pyfaults.simXRD as xs
    
    exptTT, exptInts = importSim(path, fn)
    
    truncTT = []
    truncInts = []
    for i in range(len(exptTT)):
        if exptTT[i] <= maxTT:
            truncTT.append(exptTT[i])
            truncInts.append(exptInts[i])
    
    truncTT = np.array(truncTT)
    truncInts = np.array(truncInts)
    exptQ = tt_to_q(truncTT, wl)
        
    exptIntsMin = truncInts - np.min(truncInts)
    exptNorm = xs.norm(exptIntsMin)
    
    return exptQ, exptNorm
        

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
        for i in range(2):
            txtSub = re.sub(strVars[i], str(sVecList[s][i]), sVecStr)
            sVecStr = txtSub
        sVecStrList.append(sVecStr)
        
    return probStrList, sVecStrList


#----------------------------------------------------------------------------------------
def genUnitCell(cellName, path, csvFile, lattice, layerList, stackDir, 
                childLayerList=None, childGenVecList=None):
    from pyfaults import importCSV
    import pyfaults.layer as lyr
    from pyfaults.unitcell import Unitcell
    
    csv = importCSV(path, csvFile)
    
    layers = lyr.getLayers(csv, lattice, layerList, stackDir)
    
    if childLayerList is not None:
        for i in childLayerList:
            name = childLayerList[i]
            parent = childGenVecList[i][0]
            vec = childGenVecList[i][1]
            
            for j in layers:
                if j.layerName == parent:
                    childLayer = lyr.j.genChildLayer(name, vec)
                    layers.append(childLayer)
    
    unitcell = Unitcell(cellName, layers, lattice)
    
    return unitcell


#----------------------------------------------------------------------------------------
def genSupercells(unitcell, nStacks, fltLayer, probList, sVecList, path):
    from pyfaults.supercell import Supercell
    from pyfaults import toCif

    cellList = []
    
    probStrList, sVecStrList = formatStr(probList, sVecList)
    
    df = pd.DataFrame()
    
    UF = Supercell(unitcell, nStacks)
    cellList.append([UF, "Unfaulted"])
    
    modelCol = ["Unfaulted"]
    sVecTxtCol = [np.nan]
    probTxtCol = [np.nan]
    sxCol = [np.nan]
    syCol = [np.nan]
    szCol = [np.nan]
    probCol = [np.nan]
    
    for p in range(len(probList)):
        for s in range(len(sVecList)):
            FLT = Supercell(unitcell, nStacks, fltLayer=fltLayer, 
                            stackVec=sVecList[s], stackProb=probList[p])
            
            cellTag = "S" + str(s+1) + "_P" + str(int(probList[p]*100))
            cellList.append([FLT, cellTag])
            
            modelCol.append(cellTag)
            sVecTxtCol.append(sVecStrList[s])
            probTxtCol.append(probStrList[p])
            sxCol.append(sVecList[s][0])
            syCol.append(sVecList[s][1])
            szCol.append(sVecList[s][2])
            probCol.append(probList[p])
            
    for c in range(len(cellList)):
        toCif((cellList[c][0]), path, cellList[c][1])
        
    df["Model"] = modelCol
    df["Stacking Vector"] = sVecTxtCol
    df["Stacking Probability"] = probTxtCol
    df["S_x"] = sxCol
    df["S_y"] = syCol
    df["S_z"] = szCol
    df["P"] = probCol
            
    return df


#----------------------------------------------------------------------------------------
def calcSims(df, path, wl, maxTT, pw):
    import pyfaults.simXRD as xs
    if os.path.exists(path + "sims/") == False:
        os.mkdir(path + "sims/")
    
    simDF = cp.deepcopy(df)
    
    simQList = []
    simDiffList = []

    for i in simDF.index:
        name = simDF["Model"][i]
        Q, ints = xs.fullSim(path, name, wl, maxTT, pw=pw)
        
        simQList.append(Q)
        simDiffList.append(xs.norm(ints))
        
        with open(path + "sims/" + name + "_sim.txt", "w") as f:
            for (Q, ints) in zip(Q, ints):
                f.write("{0} {1}\n".format(Q, ints))
        f.close()
        
    simDF["Simulated Q"] = simQList
    simDF["Simulated Intensity"] = simDiffList
    
    return simDF


#----------------------------------------------------------------------------------------
def calcDiffs(path, simDF, expt):
    import pyfaults.simXRD as xs
    
    if os.path.exists(path + "diffCurves/") == False:
        os.mkdir(path + "diffCurves/")
    
    exptDiffDF = cp.deepcopy(simDF)
    exptQ = expt[0]
    exptInts = expt[1]
    
    exptDiffQ = []
    exptDiffInts = []
    for row in simDF.index:
        name = simDF["Model"][row]
        simQ = simDF["Simulated Q"][row]
        simInts = simDF["Simulated Intensity"][row]
        
        diffQ, diffInts = xs.diffCurve(exptQ, simQ, exptInts, simInts)
        
        exptDiffQ.append(diffQ)
        exptDiffInts.append(diffInts)
        
        with open(path + "diffCurves/" + name + "_exptDiff.txt", "w") as f:
            for (q, ints) in zip(diffQ, diffInts):
                f.write("{0} {1}\n".format(q, ints))
        f.close() 
    
    exptDiffDF["Expt vs. Model Q"] = exptDiffQ
    exptDiffDF["Expt vs. Model Difference"] = exptDiffInts
    
    return exptDiffDF


#----------------------------------------------------------------------------------------
def calcFitDiffs(path, exptDiffDF):
    if os.path.exists(path + "fitDiffCurves/") == False:
        os.mkdir(path + "fitDiffCurves/")
    
    fitDiffDF = cp.deepcopy(exptDiffDF)
    
    UFdiff = exptDiffDF["Expt vs. Model Difference"][0]
    
    fitDiffs = []
    for i in exptDiffDF.index:
        if i == 0:
            fitDiffs.append(np.nan)
        elif i > 0:
            name = exptDiffDF["Model"][i]
            modelDiff = exptDiffDF["Expt vs. Model Difference"][i]
            
            diff = np.subtract(UFdiff, modelDiff)
                
            fitDiffs.append(diff)
            
            with open(path + "fitDiffCurves/" + name + "_fitDiff.txt", "w") as f:
                for diff in range(len(fitDiffs)):
                    f.write("{0}\n".format(diff))
            f.close() 
            
    fitDiffDF["UF vs. FLT Model"] = fitDiffs
    
    return fitDiffDF


#----------------------------------------------------------------------------------------
def makePeakDF(labels, qVals, pw):
    df = pd.DataFrame()
    
    peakRanges = []
    for i in range(len(labels)):
        qMin = float("%.3f"%(qVals[i] - (pw/2)))
        qMax = float("%.3f"%(qVals[i] + (pw/2)))
        qRange = [qMin, qMax]
        peakRanges.append(qRange)
        
    df["Reflection"] = labels
    df["Q"] = qVals
    df["Peak Range"] = peakRanges
    
    return df


#----------------------------------------------------------------------------------------
def peakFitCompare(fitDiffDF, peakDF):
    
    Q = fitDiffDF["Simulated Q"][0]
    
    fitComp = []
    
    # for each model
    for i in fitDiffDF.index:
        name = fitDiffDF["Model"][i]
        
        compDF = pd.DataFrame()
        
        Q = fitDiffDF["Expt vs. Model Q"][i]
        fullFitDiff = fitDiffDF["UF vs. FLT Model"][i]
        
        truncFitCol = []
        sumCol = []
        resultCol = []
        
        if i == 0:
            truncFitCol.append(np.nan)
            sumCol.append(np.nan)
            resultCol.append(np.nan)
        
        if i > 0:
            Q = fitDiffDF["Expt vs. Model Q"][i]
            fullFitDiff = fitDiffDF["UF vs. FLT Model"][i]
            
            # for each peak
            for j in peakDF.index:
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
                
                truncFitCol.append(truncFitDiff)
                sumCol.append(diffSum)
                resultCol.append(result)
            
            compDF["Peak"] = peakDF["Reflection"]
            compDF["UF vs. FLT Fit Difference"] = truncFitCol
            compDF["Sum"] = sumCol
            compDF["Fit Result"] = resultCol
            
            addEntry = [name, compDF]
            fitComp.append(addEntry)
        
    return fitComp


#----------------------------------------------------------------------------------------
def displayResults(fitComp):
    for i in range(len(fitComp)):
        modelResult = fitComp[i][1]
        
        print("Model: " + fitComp[i][0])
    
        for j in modelResult.index:
            print("Peak: " + modelResult["Peak"][j] + " " +\
                  modelResult["Fit Result"][j])
                
    return


#----------------------------------------------------------------------------------------
def plotSearchResults(fitDiffDF, peakDF, xSpacing, wl, yLim, rowAdj=0.0,
                      colAdj=0.0):
    from pyfaults.plotXRD import gradientGen2D

    rows = len(fitDiffDF)-1
    cols = len(peakDF)
    
    g = gradientGen2D(["#00C6BF", "#009AE1", "#5D7AD3", "#B430C2"], rows, cols)
    
    fig, (p) = plt.subplots(rows, cols, figsize=(rows*2, cols))
    
    Q = fitDiffDF["Expt vs. Model Q"][0]
    
    for r in range(rows):
        for c in range(cols):
            ints = fitDiffDF["UF vs. FLT Model"][r+1]
            p[r][c].plot(Q, ints, color=g[r][c].hex)
        
    xLims = []
    for i in peakDF.index:
        qRange = peakDF["Peak Range"][i]
        qMin = qRange[0] - xSpacing
        qMax = qRange[1] + xSpacing
        xLims.append([qMin, qMax])
        
    for r in range(rows):
        for c in range(cols):
            p[r][c].set_xlim(xLims[c][0], xLims[c][1])
            p[r][c].set_ylim(yLim)
            p[r][c].tick_params(axis="both", labelsize="14")
            
            if r != (rows-1):
                p[r][c].get_xaxis().set_visible(False)
            if c != 0:
                p[r][c].get_yaxis().set_visible(False)
            
    xLabel = r"Q (\AA" r"$^{-1}$, $\lambda=$" + str(wl) + r" \AA)"
    yLabel = r"Diff$_{\mathrm{UF}} -$ Diff$_{\mathrm{F}}$ (counts, normalized)"
    
    fig.supxlabel(xLabel, fontsize=16)
    fig.supylabel(yLabel, fontsize=16)
    
    yMid = ((yLim[1] - yLim[0]) / 2) + yLim[0]
    xEnd = xLims[-1][1]
    
    rowLabels = []
    for r in range(rows):
        sVec = str(fitDiffDF["Stacking Vector"][r+1])
        prob = str(fitDiffDF["Stacking Probability"][r+1])
        label = "\n".join((sVec, prob))
        rowLabels.append(label)
    
    colLabels = peakDF["Reflection"]
    
    for r in range(rows):
        p[r][-1].text(xEnd+rowAdj, yMid, rowLabels[r], color=g[r][-1].hex, 
                      fontsize="14", ha="left", va="center")
        
    for c in range(cols):
        xMid = ((xLims[c][1] - xLims[c][0]) / 2) + xLims[c][0]
        p[0][c].text(xMid, yLim[1]+colAdj, colLabels[c], color=g[0][c].hex, 
                     fontsize="14", ha="center", va="bottom")
        
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    
    return(p)


#----------------------------------------------------------------------------------------
def autoSearch(path, unitcell, expt, nStacks, fltLayer, probList, sVecList, 
               peakLabels, peakQ, calcPW, wl, maxTT, simPW=0.0, showPlot=False):
    
    # generate supercells
    cellDF = genSupercells(unitcell, nStacks, fltLayer, probList, sVecList, path)
    print("Finished generating supercell CIFs: ")
    for i in cellDF.index:
        print(cellDF["Model"][i])
    
    # simulate XRD
    simDF = calcSims(cellDF, path, wl, maxTT, simPW)
    print("Finished simulating X-ray diffraction patterns")
    
    # calculate expt vs. model difference curves
    exptDiffDF = calcDiffs(path, simDF, expt)
    print("Finished calculating experiment vs. model difference curves")
    
    # calculate fit differences
    fitDiffDF = calcFitDiffs(path, exptDiffDF)
    print(r"Finished calculating fit difference curves,\
          Diff$_\mathrm{UF} -$ Diff$_\mathrm{FLT}$")
    
    # set peak parameters
    peakDF = makePeakDF(peakLabels, peakQ, calcPW)
    
    # compare peak fits
    fitComp = peakFitCompare(fitDiffDF, peakDF)
    
    displayResults(fitComp)
    
    return fitComp, fitDiffDF, peakDF
















        
        