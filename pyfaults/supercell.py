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
    conType (str) : type of supercell construction, can be 'Displacement' or 'Intercalation'
    fltLayer (str, optional) : faulted layer name
    stackVec (array_like, optional) : displacement vector [x,y,z] for faulted layer in fractional coordinates
    stackProb (float, optional) : stacking fault probability
    zAdj (float, optional) :
    intLayer (Layer, optional) : 
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
                 doc='int : number of stacks in supercell')
    
    conType =\
        property(lambda self: self._conType,
                 doc='str : Displacement, Transition Matrix, or Intercalation')
    
    layers =\
        property(lambda self: self._layers, 
                 doc='list (Layer) : layers contained in supercell')
    
    fltLayer =\
        property(lambda self: self._fltLayer,
                 doc='fltLayer [optional] : faulted layer name')
    
    stackVec =\
        property(lambda self: self._stackVec, 
                 doc='nparray [optional] : stacking vector')
    
    stackProb =\
        property(lambda self: self._stackProb, 
                 lambda self, val: self.setParam(stackProb=val),
                 doc='float [optional] : stacking probability')
        
    zAdj =\
        property(lambda self: self._zAdj, 
                 lambda self, val: self.setParam(zAdj=val),
                 doc='float [optional] : z-value adjustment')
        
    intLayer =\
        property(lambda self: self._intLayer, 
                 lambda self, val: self.setLayers(intLayer=val),
                 doc='list (Layer) : inserted intercalation layer')
    
    # initialization, defines Supercell defaults ----------
    def __init__(self, unitcell, nStacks, conType, fltLayer=None, stackVec=None, stackProb=None, zAdj=None, intLayer=None):
        
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
        self._conType = None
        self._layers = None
        self._fltLayer = None
        self._stackVec = None
        self._stackProb = None
        self._zAdj = None
        self._intLayer = None
        self.setParam(nStacks, conType)
        
        if self.conType == 'Displacement':
            self.setDisLayers(unitcell, fltLayer, stackVec, stackProb)
            
        if self.conType == 'Intercalation':
            self.setIntLayers(unitcell, fltLayer, stackProb, zAdj, intLayer)
        return
    
    # sets number of stacks and construction type ----------
    def setParam(self, nStacks=None, conType=None):
        if nStacks is not None:
            self._nStacks = nStacks
        if conType is not None:
            self._conType = conType
        return
    
    # constructs displacement supercell layers ----------
    def setDisLayers(self, unitcell, fltLayer=None, stackVec=None, stackProb=None):
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

    # constructs intercalation supercell layers ----------
    def setIntLayers(self, unitcell, fltLayer=None, stackProb=None, zAdj=None, intLayer=None):
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
