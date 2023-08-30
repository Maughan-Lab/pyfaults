#########################################################################################
# pyfaults.childLayer
# Author: Sinclair R. Combs
#########################################################################################

import copy as cp
import numpy as np
from pyfaults.layerAtom import LayerAtom

''' Child object class -- Generates translated child layer from parent layer'''
#----------------------------------------------------------------------------------------

class ChildLayer (object):
    '''
    Parameters
    ----------
    layerName : str
        Unique identifier for Child object
    parent : Layer
        Generator parent layer
    transVec : nparray
        Translation vector to generate child from parent
    '''
    
    # properties ------------------------------------------------------------------------
    layerName = property(lambda self: self._layerName,
                         lambda self, val: self.setParam(layerName=val))
    
    parent = property(lambda self: self._parent,
                      lambda self, val: self.setParam(parent=val))
    
    transVec = property(lambda self: self._transVec,
                        lambda self, val: self.setParam(transVec=val))
    
    lattice = property(lambda self: self._lattice,
                       lambda self, val: self.setParam(lattice=val))
    
    # creates instance of Child object --------------------------------------------------
    def __init__(self, layerName, parent, transVec):
        
        # initialize parameters
        self._layerName = None
        self._parent = None
        self._transVec = None
        
        self._atoms = []
        self._lattice = None
        
        self.setLayerName(layerName)
        
        self.generate(parent, transVec)
        
        return
    
    # sets layer parameters -------------------------------------------------------------
    def setParam(self, layerName=None, parent=None, transVec=None):
        if layerName is not None:
            self._layerName = layerName
        if parent is not None:
            self._parent = parent
            latt = parent.lattice
            self._lattice = latt
        if transVec is not None:
            self._transVec = transVec
         
        return
    
    # creates child layer by translating parent layer -----------------------------------
    def generate(self, parent, transVec):
        alist = parent.atoms
        for i in range(len(alist)):
            parentAtom = cp.deepcopy(alist[i])
            
            splitLabel = parentAtom.atomLabel.split("_")
            newLabel = splitLabel[0] + "_" + self.layerName
            
            newPos = np.add(parentAtom.xyz, transVec)
            
            childAtom = LayerAtom(parentAtom.layerName, newLabel,
                                  parentAtom.element, newPos, 
                                  parentAtom.occupancy, parentAtom.lattice)
            
            self._atoms.append(childAtom)
            return
    
    # prints layer name and atoms -------------------------------------------------------
    def display(self):
        lyr = "\n".join(("Layer Name: " + self.layerName, 
                       "Atoms: " + str(self.atoms)))
        print(lyr)
    
    # prints atom labels and positions ---------------------------------------------------
    def atom_info(self):
        alist = self.atoms
        for i in range(len(alist)):
            print(alist[i].atomLabel + ", " + str(alist[i].xyz))