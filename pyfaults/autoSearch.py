#########################################################################################
# pyfaults.autoSearch
# Author: Sinclair R. Combs
#########################################################################################

import pandas as pd

class AutoSearch(object):
    
    # properties ------------------------------------------------------------------------
    expt = property(lambda self: self._expt)
    
    unitcell = property(lambda self: self._unitcell,
                        lambda self, val: self.setUnitCell(unitcell=val))
    
    nStacks = property(lambda self: self._nStacks,
                       lambda self, val: self.setSearchParams(nStacks=val))
    
    fltLayer = property(lambda self: self._fltLayer,
                        lambda self, val: self.setSearchParams(fltLayer=val))
    
    sVecList = property(lambda self: self._sVecList,
                        lambda self, val: self.setSearchParams(sVecList=val))
    
    sProbList = property(lambda self: self._sProbList,
                         lambda self, val: self.setSearchParams(sProbList=val))
    
    maxQ = property(lambda self: self._maxQ,
                    lambda self, val: self.setSearchParams(maxQ=val))
    
    wl = property(lambda self: self._wl,
                  lambda self, val: self.setSearchParams(wl=val))
    
    saveDir = property(lambda self: self._saveDir,
                       lambda self, val: self.setSearchParams(saveDir=val))
    
    cellList = property(lambda self: self._cellList,
                        lambda self, val: self.genSupercells(cellList=val))
    
    simList = property(lambda self: self._simList,
                       lambda self, val: self.calcSims(simList=val))
    
    exptDiffList = property(lambda self: self._exptDiffList,
                            lambda self, val: self.calcDiffs(exptDiffList=val))
    
    fitDiffList = property(lambda self: self._fitDiffList,
                            lambda self, val: self.fitDiffs(fitDiffList=val))
    
    def __init__(self, expt):
        
        # initialize parameters
        self._expt = expt
        self._unitcell = None
        self._nStacks = None
        self._fltLayer = None
        self._sVecList = None
        self._sProbList = None
        self._maxQ = None
        self._wl = None
        self._saveDir = None
        
        self._cellList = None
        self._cellList = self.genSupercells()
        
        self._simList = None
        self._simList = self.calcSims()
        
        self._exptDiffList = None
        self._exptDiffList = self.calcDiffs()
        
        self._fitDiffList = None
        self._fitDiffList = self.fitDiffs()
        
        return
    
    # sets unit cell from CSV file ------------------------------------------------------
    def setUnitCell(self, path, fn, lattice, layerNames, stackDir):
        from pyfaults.lattice import Lattice
        from pyfaults.layerAtom import LayerAtom
        from pyfaults.layer import Layer
        from pyfaults.unitcell import Unitcell
        
        df = pd.read_csv(path + fn + ".csv")
        
        # set lattice parameters so stacking direction is along c
        if stackDir == "a":
            newLatt = Lattice(a=lattice.c, b=lattice.b, c=lattice.a,
                              alpha=lattice.gamma, beta=lattice.beta, 
                              gamma=lattice.alpha)
        elif stackDir == "b":
            newLatt = Lattice(a=lattice.a, b=lattice.c, c=lattice.b,
                              alpha=lattice.alpha, beta=lattice.gamma, 
                              gamma=lattice.beta)
        elif stackDir == "c":
            newLatt = Lattice(a=lattice.a, b=lattice.b, c=lattice.c,
                              alpha=lattice.alpha, beta=lattice.beta, 
                              gamma=lattice.gamma)
        
        # create layer objects
        layers = []
        for i in layerNames:
            lyr = df[df["Layer"] == i]
            
            alist = []
            for j in range(len(lyr.index)):
                layerName = df.iloc[j]["Layer"]
                atomLabel = df.iloc[j]["Atom"]
                element = df.iloc[j]["Element"]
                x = df.iloc[j]["x"]
                y = df.iloc[j]["y"]
                z = df.iloc[j]["z"]
                occ = df.iloc[j]["Occupancy"]
                
                if stackDir == "a":
                    newXYZ = [z, y, x]
                elif stackDir == "b":
                    newXYZ = [x, z, y]
                elif stackDir == "c":
                    newXYZ = [x, y, z]
                
                newAtom = LayerAtom(layerName, atomLabel, element, newXYZ, occ, newLatt)
                alist.append(newAtom)    
        
            newLayer = Layer(alist, newLatt, i)
            layers.append(newLayer)
        
        # create unit cell object
        unitcell = Unitcell(fn, layers, newLatt)
        
        self._unitcell = unitcell
        
        return

    # sets search parameters ------------------------------------------------------------
    def setSearchParams(self, nStacks, fltLayer, sVecList, sProbList, maxQ, wl, saveDir):
        
        self._nStacks = nStacks
        self._fltLayer = fltLayer
        self._sVecList = sVecList
        self._sProbList = sProbList
        self._maxQ = maxQ
        self._wl = wl
        self._saveDir = saveDir
        
        return
    
    # generate supercells ---------------------------------------------------------------
    def genSupercells(self):
        from pyfaults import toCif
        from pyfaults.supercell import Supercell
        
        if self.nStacks is None:
            raise ValueError("Search parameters not specified")
            
        cellList = []
            
        UFcell = Supercell(self.unitcell, self.nStacks)
        cellList.append([UFcell, "Unfaulted", 0, "UF"])
        toCif(UFcell, self.saveDir, "UF")
        
        for i in range(len(self.sProbList)):
            for j in range(len(self.sVecList)):
                FLTcell = Supercell(self.unitell, self.nStacks, 
                                    fltLayer=self.fltLayer, stackVec=self.sVecList[j], 
                                    stackProb=self.sProbList[i])
                
                sx = str(self.sVecList[j][0])
                sy = str(self.sVecList[j][1])
                sz = str(self.sVecList[j][2])
                
                sVec = "[" + sx + ", " + sy + ", " + sz + "]"
                    
                tag = "P" + str(i) + "_FLT_" + str(j)
                
                cellList.append([FLTcell, sVec, self.sProbList[i], tag])
                
                toCif(FLTcell, self.saveDir, tag)
                
        self._cellList = cellList
        
        return
    
    # calculate XRD simulations ---------------------------------------------------------
    def calcSims(self):
        from pyfaults import q_to_tt
        import pyfaults.simXRD as xs
        
        if self.cellList is None:
            raise ValueError("Supercells not specified")
            
        ttMax = q_to_tt(self.maxQ, self.wl)
            
        simList = []
        
        for i in range(len(self.cellList)):
            cif = self.cellList[i][3]
            
            Q, ints = xs.fullSim(self.saveDir, cif, self.wl, ttMax, pw=0.01)
            
            simList.append([Q, ints])
            
        self._simList = simList
        
        return
    
    # calculate model vs experimental differences ---------------------------------------
    def calcDiffs(self):
        from pyfaults.simXRD import diffCurve
        
        if self.simList is None:
            raise ValueError("Diffraction simulations not calculated")
        
        exptDiffList = []
        
        for i in range(len(self.simList)):
            diffQ, diffInts = diffCurve(self.expt[0], self.simList[i][0], 
                                        self.expt[1], self.simList[i][1])
            exptDiffList.append([diffQ, diffInts])
            
        self._exptDiffList = exptDiffList
        
        return
    
    # calculate goodness of fits --------------------------------------------------------
    def fitDiffs(self):
        from pyfaults.simXRD import diffCurve
        
        if self.exptDiffList is None:
            raise ValueError("Difference curves not calculated")
        
        fitDiffList = []
        
        for i in range(1, len(self.exptDiffList)):
            diffQ, diffInts = diffCurve(self.exptDiffList[0][0], self.exptDiffList[i][0],
                                        self.exptDiffList[0][1], self.exptDiffList[i][1])
            
        self._fitDiffList = fitDiffList
        
        return


















        
        