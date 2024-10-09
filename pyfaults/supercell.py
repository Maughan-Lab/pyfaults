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
    def __init__(self, unitcell, nStacks, fltLayer=None, stackVec=None, stackProb=None, zAdj=None, intLayer=None):
        
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
        self._zAdj = None
        self._intLayer = None
        self.setParam(nStacks)
        
        self.setLayers(unitcell, fltLayer, stackVec, stackProb, zAdj, intLayer)
        return
    
    # sets number of stacks and construction type ----------
    def setParam(self, nStacks=None):
        if nStacks is not None:
            self._nStacks = nStacks
        return
    
    # constructs supercell layers ----------
    def setLayers(self, unitcell, fltLayer=None, stackVec=None, stackProb=None, zAdj=None, intLayer=None):
        if fltLayer is not None:
            self._fltLayer = fltLayer
        if stackVec is not None:
            self._stackVec = stackVec
        if stackProb is not None:
            self._stackProb = stackProb
        if zAdj is not None:
            self._zAdj = zAdj
        if intLayer is not None:
            self._intLayer = intLayer
        
        p = []
        for i in range(self.nStacks):
            p.append(r.randint(0,100))
            
        fltCount = 0
        for i in range(len(p)):
            if p[i] <= (stackProb * 100):
                fltCount = fltCount + 1
                
        if zAdj is not None:
            self.lattice.c = self.lattice.c + (zAdj*fltCount)
        
        newLayers = []
        for n in range(self.nStacks):
            tag = '_n' + str(n+1)
            for lyr in self.unitcell.layers:
                    
                if lyr.layerName == self.fltLayer:
                    if p[n] <= (stackProb * 100):
                        newFltLyr = cp.deepcopy(lyr)
                        newLayerName = lyr.layerName + tag + '_fault'
                        
                        if zAdj is None:
                            newFltLyr.setParam(layerName=newLayerName, lattice=self.lattice)
                            
                            for atom in newFltLyr.atoms:
                                alabel = atom.atomLabel.split('_')
                                newXYZ = [atom.x, atom.y, ((atom.z + n) / self.nStacks)]
                                fltXYZ = np.add(newXYZ, stackVec)
                                atom.setParam(layerName=newLayerName, atomLabel=alabel[0], xyz=fltXYZ, lattice=self.lattice)
                            
                            newLayers.append(newFltLyr)
                            
                        if zAdj is not None:
                            self.lattice.c = self.lattice.c + zAdj
                                
                            newFltLyr.setParam(layerName=newLayerName, lattice=self.lattice)
                            
                            for atom in newFltLyr.atoms:
                                alabel = atom.atomLabel.split('_')
                                newXYZ = [atom.x, atom.y, ((atom.z + n) / self.nStacks)]
                                fltXYZ = np.add(newXYZ, [0,0,zAdj])
                                atom.setParam(layerName=newLayerName, atomLabel=alabel[0], xyz=fltXYZ, lattice=self.lattice)
                            
                            newLayers.append(newFltLyr)

                            if self.intLayer is not None:
                                newIntLayer = self.addIntLayer(n, tag)
                                newLayers.append(newIntLayer)
                        
                    
                    else:
                        newUFLyr = cp.deepcopy(lyr)
                        
                        newLayerName = lyr.layerName + tag
                        newUFLyr.setParam(layerName=newLayerName, lattice=self.lattice)
                        for atom in newUFLyr.atoms:
                            alabel = atom.atomLabel.split('_')
                            newXYZ = [atom.x, atom.y, ((atom.z + n) / self.nStacks)]
                            atom.setParam(layerName=newLayerName, atomLabel=alabel[0], xyz=newXYZ, lattice=self.lattice)
                
                    newLayers.append(newUFLyr)
                    
        self._layers = newLayers
        return
    
    # adds intercalation layer ----------
    def addIntLayer(self, n, tag):
        newLyr = cp.deepcopy(self.intLayer)
        newLyr.setParam(layerName='I' + tag, lattice=self.lattice)
        
        for atom in newLyr.atoms:
            alabel = atom.atomLabel.split('_')
            newXYZ = [atom.x, atom.y, ((atom.z + n) / self.nStacks)]
            
            atom.setParam(layerName='I' + tag, atomLabel=alabel[0], xyz=newXYZ, lattice=self.lattice)
            
        return newLyr
    
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
