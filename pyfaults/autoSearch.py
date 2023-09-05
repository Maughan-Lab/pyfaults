#########################################################################################
# pyfaults.autoSearch
# Author: Sinclair R. Combs
#########################################################################################

class AutoSearch(object):
    
    # properties ------------------------------------------------------------------------
    expt = property(lambda self: self._expt)
    
    unitcell = property(lambda self: self._unitcell)
    
    nStacks = property(lambda self: self._nStacks)
    
    fltLayer = property(lambda self: self._fltLayer)
    
    sVecList = property(lambda self: self._sVecList)
    
    sProbList = property(lambda self: self._sProbList)
    
    maxQ = property(lambda self: self._maxQ)
    
    wl = property(lambda self: self._wl)
    
    saveDir = property(lambda self: self._saveDir)
    
    cellList = property(lambda self: self._cellList,
                        lambda self, val: self.genSupercells(cellList=val))
    
    simList = property(lambda self: self._simList,
                       lambda self, val: self.calcSims(simList=val))
    
    exptDiffList = property(lambda self: self._exptDiffList,
                            lambda self, val: self.calcDiffs(exptDiffList=val))
    
    fitDiffList = property(lambda self: self._fitDiffList,
                            lambda self, val: self.fitDiffs(fitDiffList=val))
    
    def __init__(self, expt, unitcell, nStacks, fltLayer, sVecList, sProbList,
                 maxQ, wl, saveDir):
        
        # initialize parameters
        self._expt = expt
        self._unitcell = unitcell
        self._nStacks = nStacks
        self._fltLayer = fltLayer
        self._sVecList = sVecList
        self._sProbList = sProbList
        self._maxQ = maxQ
        self._wl = wl
        self._saveDir = saveDir
        
        self._cellList = self.genSupercells()
        
        self._simList = self.calcSims()
        
        self._exptDiffList = self.calcDiffs()
        
        self._fitDiffList = self.fitDiffs()
        
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
                FLTcell = Supercell(self.unitcell, self.nStacks, 
                                    fltLayer=self.fltLayer, stackVec=self.sVecList[j], 
                                    stackProb=self.sProbList[i])
                
                sx = str(self.sVecList[j][0])
                sy = str(self.sVecList[j][1])
                sz = str(self.sVecList[j][2])
                
                sVec = "[" + sx + ", " + sy + ", " + sz + "]"
                    
                tag = "P" + str(i) + "_FLT_" + str(j)
                
                cellList.append([FLTcell, sVec, self.sProbList[i], tag])
                
                toCif(FLTcell, self.saveDir, tag)
        
        return cellList
    
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
        
        return simList
    
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
        
        return exptDiffList
    
    # calculate goodness of fits --------------------------------------------------------
    def fitDiffs(self):
        from pyfaults.simXRD import diffCurve
        
        if self.exptDiffList is None:
            raise ValueError("Difference curves not calculated")
        
        fitDiffList = []
        
        for i in range(1, len(self.exptDiffList)):
            diffQ, diffInts = diffCurve(self.exptDiffList[0][0], self.exptDiffList[i][0],
                                        self.exptDiffList[0][1], self.exptDiffList[i][1])
        
        return fitDiffList


















        
        