#########################################################################################
# pyfaults.supercell
# Author: Sinclair R. Combs
#########################################################################################

import copy as cp
import numpy as np
import random as r

''' Supercell object class -- Stores supercell properties '''
#----------------------------------------------------------------------------------------

class Supercell(object):
    '''
    Parameters
    ----------
    unitcell : Unitcell
        Base unit cell
    nStacks : int
        Number of stacks / unit cells in supercell
    fltLayer : str, optional
        Name of faulted layer. The default is None.
    stackVec : nparray, optional
        Stacking vector relative to ideal position. The default is None.
    stackProb : float, optional
        Stacking probability. The default is None.
    '''
    
    # properties ------------------------------------------------------------------------
    unitcell = property(lambda self: self._unitcell)
    
    lattice = property(lambda self: self._lattice)
    
    nStacks = property(lambda self: self._nStacks, 
                       lambda self, val: self.setParam(nStacks=val))
    
    layers = property(lambda self: self._layers, 
                      lambda self, val: self.setLayers(layers=val))
    
    fltLayer = property(lambda self: self._fltLayer,
                        lambda self, val: self.setParam(fltLayer=val))
    
    stackVec = property(lambda self: self._stackVec, 
                    lambda self, val: self.setParam(stackVec=val))
    
    stackProb = property(lambda self: self._stackProb, 
                    lambda self, val: self.setParam(stackProb=val))
    
    # creates instance of Supercell object ----------------------------------------------
    def __init__(self, unitcell, nStacks, fltLayer=None, stackVec=None, stackProb=None):
        from pyfaults.lattice import Lattice

        # initialize parameters
        self._unitcell = unitcell
        
        newLatt = Lattice(unitcell.lattice.a,
                          unitcell.lattice.b,
                          (unitcell.lattice.c * nStacks),
                          unitcell.lattice.alpha,
                          unitcell.lattice.beta,
                          unitcell.lattice.gamma)
        
        self._lattice = newLatt
        
        self._nStacks = None
        self._layers = None
        
        self._fltLayer = None
        self._stackVec = None
        self._stackProb = None
        
        self.setParam(nStacks)
        self.setLayers(unitcell, fltLayer, stackVec, stackProb)
        
        return
    
    # set cell parameters ---------------------------------------------------------------
    def setParam(self, nStacks=None):
        
        if nStacks is not None:
            self._nStacks = nStacks
            
        return
    
    # define layers of supercell --------------------------------------------------------
    def setLayers(self, unitcell, fltLayer=None, stackVec=None, stackProb=None):
        
        newLayers = []
        
        for lyr in self.unitcell.layers:
            
            for n in range(self.nStacks):
                
                tag = "_n" + str(n+1)
                p = r.randint(0,100)
                
                if fltLayer is not None:
                    if lyr.layerName == fltLayer:
                        if p <= (stackProb * 100):
                            isFlt = True
                        elif p > (stackProb * 100):
                            isFlt = False
                    elif lyr.layerName != fltLayer:
                        isFlt = False
                elif fltLayer is None: 
                    isFlt = False
                
                newLyr = cp.deepcopy(lyr)
                
                if isFlt == True:
                    newLayerName = lyr.layerName + tag + "_fault"
                    newLyr.setParam(layerName=newLayerName, lattice=self.lattice)
                    
                    for atom in newLyr.atoms:
                        alabel = atom.atomLabel.split("_")
                        
                        newXYZ = [atom.x, atom.y, ((atom.z + n) / self.nStacks)]
                        fltXYZ = np.add(newXYZ, stackVec)
                        
                        atom.setParam(layerName=newLayerName, 
                                      atomLabel=alabel[0], xyz=fltXYZ, 
                                      lattice=self.lattice)
                if isFlt == False:
                    newLayerName = lyr.layerName + tag
                    newLyr.setParam(layerName=newLayerName, lattice=self.lattice)
                    
                    for atom in newLyr.atoms:
                        alabel = atom.atomLabel.split("_")
                        
                        newXYZ = [atom.x, atom.y, ((atom.z + n) / self.nStacks)]
                        
                        atom.setParam(layerName=newLayerName, 
                                      atomLabel=alabel[0], xyz=newXYZ, 
                                      lattice=self.lattice)
                
                newLayers.append(newLyr)
        
                
        self._layers = newLayers
        
        return
    
    # prints layer names ----------------------------------------------------------------
    def layer_info(self):
        for i in self.layers:
            print(i.layerName)
            
    def show_faults(self):
        for lyr in self.layers:
            if "fault" in lyr.layerName:
                print(lyr.layerName) 