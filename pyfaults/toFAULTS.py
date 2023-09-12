#########################################################################################
# pyfaults.toFAULTS
# Author: Sinclair R. Combs
#########################################################################################

import pandas as pd

def toXYSigma(path, fn, Q, ints):
    with open(path + fn + ".dat", "w") as f:
        f.write("XYSIGMA \n fn \n -----------------------------------")
        for (Q, ints) in zip(Q, ints):
            f.write("{0} {1} 0.0 \n").format(Q, ints)
    f.close()
    

def toFAULTS(title, wl, instBroad, unitcell, fltLyr, sVec, fltProb, t, ttMin, 
             exptFile=None, bgCoeff=None):
    # instBroad = [type, u, v, w, x, Dg, Dl]
    # t= layer order
    
    lines = []
    
    # title section ---------------------------------------------------------------------
    lines.extend(["TITLE", title, ""])
    
    # instrumental section --------------------------------------------------------------
    if instBroad[0] == "Gaussian":
        u = str(instBroad[1])
        v = str(instBroad[2])
        w = str(instBroad[3])
        x = "0.0"
        Dg = str(instBroad[4])
        Dl = "0.0"
    elif instBroad[0] == "Lorentzian":
        u = "0.0"
        v = "0.0"
        w = "0.0"
        x = str(instBroad[1])
        Dg = "0.0"
        Dl = str(instBroad[2])
    elif instBroad[0] == "Pseudo-Voigt":
        u = str(instBroad[1])
        v = str(instBroad[2])
        w = str(instBroad[3])
        x = str(instBroad[4])
        Dg = str(instBroad[5])
        Dl = str(instBroad[6])
    
    instBroadStr = instBroad[0] + "\t" + u + "\t" + v + "\t" + w + "\t" + x +\
        "\t" + Dg + "\t" + Dl + "\t Trim"
    
    lines.extend(["Instrumental And Size Broadening",
                  "Radiation \t X-Ray",
                  "! \t\t lambda1 \t lambda2 \t ratio",
                  "Wavelength \t" + str(wl) + "\t 0.0 \t\t 0.0",
                  "! Instrumental aberrations",
                  "Aberrations \t 0.0000 \t 0.0000 \t 0.0000",
                  "\t \t 0.00 \t\t 0.00 \t\t 0.00",
                  "! Broadening \t u \t v \t w \t x \t Dg \t Dl",
                  instBroadStr,
                  "\t \t 0.00 \t 0.00 \t 0.00 \t 0.00 \t 0.00 \t 0.00",
                  ""])
    
    # structural section ----------------------------------------------------------------
    a = "%.6g" % (unitcell.lattice.a)
    b = "%.6g" % (unitcell.lattice.b)
    c = "%.6g" % (unitcell.lattice.c)
    alpha = "%.6g" % (unitcell.lattice.alpha)
    beta = "%.6g" % (unitcell.lattice.beta)
    gamma = "%.6g" % (unitcell.lattice.gamma)
    
    lattStr = a + "\t" + b + "\t" + c + "\t" + alpha + "\t" + beta + "\t" + gamma
    cellStr = a + "\t" + b + "\t" + c + "\t" + gamma
    
    lines.extend(["Structural",
                  "SPGR \t P1",
                  "! \t a \t b \t c \t alpha \t beta \t gamma",
                  "Avercell " + lattStr,
                  "! \t a \t b \t c \t gamma",
                  "Cell " + cellStr])
    
    numLayerTypes = len(unitcell.layers)
    
    lines.extend(["! Laue symmetry",
                  "Symm unknown",
                  "! Number of layer types",
                  "Nlayers " + str(numLayerTypes),
                  "! Layer width",
                  "Lwidth infinite",
                  ""])
    
    # layer section ---------------------------------------------------------------------
    df = pd.DataFrame()
    lyrList = []
    atomCol = []
    
    for i in range(len(unitcell.layers)):
        lyrList.append(unitcell.layers[i].layerName)
        
        lyrAtoms = []
        for j in range(len(unitcell.layers[i].atoms)):
            a = unitcell.layers[i].atoms[j]
            
            aVals = [a.element, str(a.x), str(a.y), str(a.z), str(a.occupancy)]
            lyrAtoms.append(aVals)
        atomCol.append(lyrAtoms)
        
    for i in range(len(unitcell.layers)):
        if unitcell.layers[i].layerName == fltLyr:
            lyrList.append(unitcell.layers[i].layerName + "_FLT")
            
            fltLyrAtoms = []
            for j in range(len(unitcell.layers[i].atoms)):
                a = unitcell.layers[i].atoms[j]
                
                newX = "%.6g" % (a.x + sVec[0])
                newY = "%.6g" % (a.y + sVec[1])
                newZ = "%.6g" % (a.z + sVec[2])
                
                fltAVals = [a.element, newX, newY, newZ, str(a.occupancy)]
                fltLyrAtoms.append(fltAVals)
                
            atomCol.append(fltLyrAtoms)
    
    df["Layer Name"] = lyrList
    df["Atoms"] = atomCol
    
    for i in df.index:
        lines.extend(["Layer " + str(i+1),
                      "! Layer symmetry",
                      "LYSM none"])
        
        for a in range(len(df["Atoms"][i])):
            atom = df["Atoms"][i][a]
            aStr = atom[0] + "\t" + str(a) + "\t" + atom[1] + "\t" + atom[2] +\
                "\t" + atom[3] + "\t" + "2.0" + "\t" + atom[4]
            
            lines.extend(["! Atom name \t num \t x \t y \t z \t Biso \t Occ",
                          "Atom " + aStr,
                          "\t \t \t 0.00 \t 0.00 \t 0.00 \t 0.00 \t 0.00"])
            
        lines.extend([""])
        
    # stacking section ------------------------------------------------------------------
    lines.extend(["Stacking",
                  "! Stacking type",
                  "Recursive",
                  "! Number of layers",
                  "Infinite",
                  ""])
    
    # transitions section ---------------------------------------------------------------
    lines.extend(["Transitions"])
    
    for i in range(len(df.index)-1):
        iMax = len(df.index)-1
        if (i+1) >= iMax:
            nextLyr = ((i+1)-iMax)
        else:
            nextLyr = (i+1)
        
        # layer i to unfaulted layers
        for j in range(len(df.index)-1):
            lines.extend(["! Layer " + str(i+1) + " to layer " + str(j+1)])
        
            if df["Layer Name"][j] == df["Layer Name"][i]:
                lines.extend(["LT \t 0 \t 0.0 \t 0.0 \t 0.0",
                              "\t 0.0 \t 0.0 \t 0.0 \t 0.0",
                              "FW \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0",
                              "\t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0"])
            
            if df["Layer Name"][j] == df["Layer Name"][nextLyr]:
                if df["Layer Name"][j] == fltLyr:
                    lines.extend(["LT \t" + str(1-fltProb) + "\t 0.0 \t 0.0 \t 0.0",
                                  "\t" + str(i+1) + str(j+1) + "\t 0.0 \t 0.0 \t 0.0"])
                else:
                    lines.extend(["LT \t 1 \t 0.0 \t 0.0 \t 0.0",
                                  "\t 0.0 \t 0.0 \t 0.0 \t 0.0"])
                    
                lines.extend(["FW \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0",
                              "\t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0"])
    
        # layer i to faulted layer
        lines.extend(["! Layer " + str(i+1) + " to layer 3"])
    
        if df["Layer Name"][i] == fltLyr:
            lines.extend(["LT \t 0.0 \t 0.0 \t 0.0 \t 0.0",
                          "\t 0.0 \t 0.0 \t 0.0 \t 0.0",
                          "FW \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0",
                          "\t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0"])
        else:
            sx = "%.6g" % (sVec[0])
            sy = "%.6g" % (sVec[1])
            sz = "%.6g" % (sVec[2])
            lines.extend(["LT \t" + str(fltProb) + "\t" + sx + "\t" + sy + "\t" + sz,
                          "\t" + str(i+1) + "3 \t 0.0 \t 0.0 \t 0.0"])
            lines.extend(["FW \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0",
                          "\t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0"])
        
    # faulted layer transitions
    nextFromFLT = ""
    for i in df.index:
        if df["Layer Name"][i] == fltLyr:
            nextFromFLT = df["Layer Name"][i+1]
        
        for j in df.index:
            lines.extend(["! Layer " + str(len(df.index)) + " to layer " + str(j+1)])
    
            if df["Layer Name"][j] != nextFromFLT:
                lines.extend(["LT \t 0.0 \t 0.0 \t 0.0 \t 0.0",
                              "\t 0.0 \t 0.0 \t 0.0 \t 0.0",
                              "FW \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0",
                              "\t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0"])
            elif df["Layer Name"][j] == nextFromFLT:
                sxn = "%.6g" % (-1 * sVec[0])
                syn = "%.6g" % (-1 * sVec[1])
                szn = "%.6g" % (-1 * sVec[2])
                lines.extend(["LT \t 1.0 \t 0.0 \t" + "\t" + sxn + "\t" + syn +\
                              "\t" + szn,
                              "\t 0.0 \t 0.0 \t 0.0 \t 0.0",
                              "FW \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0",
                              "\t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0"])
        
    lines.extend([""])
        
        
    # calculation section ---------------------------------------------------------------
    if exptFile is None:
        lines.extend(["Calculation",
                      "Simulation",
                      "Powder " + str(ttMin) + " 90 0.02",
                      "Replace_Files",
                      ""])
        
    if exptFile is not None:
        lines.extend(["Calculation",
                      "LMA",
                      "! Maximum correlation parameter",
                      "Corrmax \t 30",
                      "! Maximum number of function evaluations",
                      "Maxfun \t 2400",
                      "! Tolerance",
                      "Tol \t 0.1E-04",
                      "Nprint \t 0",
                      "Replace_Files",
                      ""])
    
    # experimental section --------------------------------------------------------------
    if exptFile is not None:
        bgCoeffStr = ""
        bgCodeStr = ""
        for i in len(bgCoeff):
            bgCoeffStr = bgCoeffStr + "\t" + str(bgCoeff[i])
            bgCodeStr = bgCodeStr + "\t 0.0"
    
        lines.extend(["Experimental",
                      "FILE " + exptFile + "\t 1.0 \t 0.00",
                      "Excluded_Regions 1",
                      "\t 0.0 \t" + str(ttMin),
                      "FFORMAT XYSIGMA",
                      "! Polynomial number of coefficients",
                      "Bgrcheb " + len(bgCoeff),
                      "! Polynomial coefficients",
                      bgCoeffStr,
                      bgCodeStr])
    
    return lines
    
    
    
    
    
    
    
    