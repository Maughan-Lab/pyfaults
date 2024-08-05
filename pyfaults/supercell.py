##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import copy as cp
import numpy as np
import random as r

# Supercell object class ----------
class Supercell(object):
'''
Parameters
----------
unitcell (Unitcell) : base unit cell
nStacks (int) : number of unit cell stacks in supercell
fltLayer (str, optional) : faulted layer name
stackVec (array_like, optional) : displacement vector [x,y,z] for faulted layer in fractional coordinates
stackProb (float, optional) : stacking fault probability
'''
    # properties ----------
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
    
    fltLayer =\
        property(lambda self: self._fltLayer,
                 lambda self, val: self.setParam(fltLayer=val),
                 doc='fltLayer [optional] : faulted layer name')
    
    stackVec =\
        property(lambda self: self._stackVec, 
                 lambda self, val: self.setParam(stackVec=val),
                 doc='nparray [optional] : stacking vector')
    
    stackProb =\
        property(lambda self: self._stackProb, 
                 lambda self, val: self.setParam(stackProb=val),
                 doc='float [optional] : stacking probability')
    
    # initialization, defines Supercell defaults ----------
    def __init__(self, unitcell, nStacks, 
                 fltLayer=None, stackVec=None, stackProb=None):
        from pyfaults.lattice import Lattice
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
    
    # sets number of stacks in Supercell ----------
    def setParam(self, nStacks=None):
        if nStacks is not None:
            self._nStacks = nStacks
        return
    
    # constructs supercells layers ----------
    def setLayers(self, unitcell, fltLayer=None, stackVec=None, stackProb=None):
        newLayers = []
        for lyr in self.unitcell.layers:
            for n in range(self.nStacks):
                tag = '_n' + str(n+1)
                p = r.randint(0,100)
                
                if fltLayer is not None:
                    self._fltLayer = fltLayer
                    self._stackVec = stackVec
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
                        fltXYZ = np.add(newXYZ, stackVec)
                        atom.setParam(layerName=newLayerName, 
                                      atomLabel=alabel[0], 
                                      xyz=fltXYZ, 
                                      lattice=self.lattice)
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
    
    # prints names of faulted layers ----------
    def show_faults(self):
        for lyr in self.layers:
            if "fault" in lyr.layerName:
                print(lyr.layerName)
        return
                
    # generates CIF of layer ----------
    def toCif(self, path):
        from pyfaults import toCif
        name = self.unitcell.name
        toCif(self, path, 'Supercell_' + name + '_N' + str(self.nStacks) + '_CIF')
