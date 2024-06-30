##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import copy as cp
import numpy as np
import random as r

#---------------------------------------------------------------------------------
# Supercell object class
#---------------------------------------------------------------------------------
class SupercellAdjZ(object):    
    # properties -----------------------------------------------------------------
    unitcell =\
        property(lambda self: self._unitcell,
                 doc='Unitcell : base unit cell')
    
    lattice =\
        property(lambda self: self._lattice,
                 doc='Lattice : unit cell lattice parameters')
    
    nStacks =\
        property(lambda self: self._nStacks, 
                 lambda self, val: self.setParam(nStacks=val),
                 doc='int : number of stacks in supercell')
    
    layers =\
        property(lambda self: self._layers, 
                 lambda self, val: self.setLayers(layers=val),
                 doc='list (Layer) : layers contained in supercell')
        
    intLayer =\
        property(lambda self: self._intLayer, 
                 lambda self, val: self.setLayers(intLayer=val),
                 doc='list (Layer) : inserted intercalation layer')
    
    fltLayer =\
        property(lambda self: self._fltLayer,
                 lambda self, val: self.setParam(fltLayer=val),
                 doc='fltLayer [optional] : faulted layer name')
    
    zAdj =\
        property(lambda self: self._zAdj, 
                 lambda self, val: self.setParam(zAdj=val),
                 doc='float [optional] : z-value adjustment')
    
    stackProb =\
        property(lambda self: self._stackProb, 
                 lambda self, val: self.setParam(stackProb=val),
                 doc='float [optional] : stacking probability')
    
    # creates instance of Supercell object ---------------------------------------
    def __init__(self, unitcell, nStacks, intLayer=None,
                 fltLayer=None, zAdj=None, stackProb=None):
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
        self._intLayer = None
        self._fltLayer = None
        self._zAdj = None
        self._stackProb = None
        self.setParam(nStacks)
        self.setLayers(unitcell, intLayer, fltLayer, zAdj, stackProb)
        return
    
    # set cell parameters --------------------------------------------------------
    def setParam(self, nStacks=None):
        if nStacks is not None:
            self._nStacks = nStacks
        return
    
    # set supercell layers --------------------------------------------------------
    def setLayers(self, unitcell, intLayer=None, fltLayer=None, zAdj=None, stackProb=None):
        newLayers = []
        for lyr in self.unitcell.layers:
            for n in range(self.nStacks):
                tag = '_n' + str(n+1)
                p = r.randint(0,100)
                
                if fltLayer is not None:
                    self._fltLayer = fltLayer
                    self._zAdj = zAdj
                    self._stackProb = stackProb
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
                    newLayerName = lyr.layerName + tag + '_fault'
                    newLyr.setParam(layerName=newLayerName, 
                                    lattice=self.lattice)
                    for atom in newLyr.atoms:
                        alabel = atom.atomLabel.split('_')
                        newXYZ = [atom.x, 
                                  atom.y, 
                                  ((atom.z + n) / self.nStacks)]
                        fltXYZ = np.add(newXYZ, [0,0,zAdj])
                        atom.setParam(layerName=newLayerName, 
                                      atomLabel=alabel[0], 
                                      xyz=fltXYZ, 
                                      lattice=self.lattice)
                        newLayers.append(intLayer)
                if isFlt == False:
                    newLayerName = lyr.layerName + tag
                    newLyr.setParam(layerName=newLayerName, 
                                    lattice=self.lattice)
                    for atom in newLyr.atoms:
                        alabel = atom.atomLabel.split('_')
                        newXYZ = [atom.x, 
                                  atom.y, 
                                  ((atom.z + n) / self.nStacks)]
                        atom.setParam(layerName=newLayerName, 
                                      atomLabel=alabel[0], 
                                      xyz=newXYZ, 
                                      lattice=self.lattice)
                
                newLayers.append(newLyr)
        self._layers = newLayers
        return
    
    # prints names of faulted layers ---------------------------------------------
    def show_faults(self):
        for lyr in self.layers:
            if "fault" in lyr.layerName:
                print(lyr.layerName)
        return
                
    # generates CIF of layer -----------------------------------------------------
    def toCif(self, path):
        from pyfaults import toCif
        name = self.unitcell.name
        toCif(self, path, 'Supercell_' + name + '_N' + str(self.nStacks) + '_CIF')