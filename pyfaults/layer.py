##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import copy as cp
import numpy as np

''' LAYER CLASS'''
#---------------------------------------------------------------------------------
class Layer(object):
    # properties -----------------------------------------------------------------
    atoms =\
        property(lambda self: self._atoms,
                 lambda self, val: self.setParam(atoms=val),
                 doc='list (LayerAtom): atoms in layer')
    
    lattice =\
        property(lambda self: self._lattice,
                 lambda self, val: self.setParam(lattice=val),
                 doc='Lattice : unit cell lattice parameters')
    
    layerName\
        = property(lambda self: self._layerName,
                   lambda self, val: self.setParam(layerName=val),
                   doc='str : unique identifier for Layer instance')
    
    # creates instance of Layer object -------------------------------------------
    def __init__(self, atoms, lattice, layerName):
        # initialize parameters
        self._layerName = None
        self._atoms = None
        self._lattice = None
        self.setParam(atoms, lattice, layerName)
        return
    
    # sets parameters ------------------------------------------------------------
    def setParam(self, atoms=None, lattice=None, layerName=None):
        if atoms is not None:
            self._atoms = atoms
        if lattice is not None:
            self._lattice = lattice
        if layerName is not None:
            self._layerName = layerName
        return
    
    # prints layer information ---------------------------------------------------
    def display(self):
        print(self.layerName)
        for i in self.atoms:
            print(i.atomLabel + ', ' + str(i.xyz))
        return 
            
    # generates CIF of layer -----------------------------------------------------
    def toCif(self, path):
        from pyfaults import toCif
        toCif(self, path, 'Layer_' + self.layerName + '_CIF')
    
    # generates a child layer ----------------------------------------------------
    def genChildLayer(self, childName, transVec):
        '''
        Parameters
        ----------
        childName
            str : unique identifier for child layer
        transVec
            nparray : translation vector [x, y, z] relative to parent

        Returns
        -------
        childLayer : Layer
        '''
        from pyfaults.layerAtom import LayerAtom
        
        childAtoms = []
        for a in self.atoms:
            # create copy of parent atom
            pAtom = cp.deepcopy(a)
            splitLabel = pAtom.atomLabel.split('_')
            
            # apply translation vector to atomic position
            newPos = np.add(pAtom.xyz, transVec)
            for i in range(len(newPos)):
                if newPos[i] >= 1:
                    newPos[i] = newPos[i] - 1
                    
            # create new LayerAtom instance
            cAtom = LayerAtom(childName, 
                              splitLabel[0], 
                              pAtom.element, 
                              newPos, 
                              pAtom.occupancy,
                              pAtom.biso,
                              self.lattice)
            childAtoms.append(cAtom)
            
        # create new Layer instance
        childLayer = Layer(childAtoms, self.lattice, childName)
        return childLayer


''' LAYER CLASS HELPER METHODS '''
#---------------------------------------------------------------------------------
# import layers from dataframe ---------------------------------------------------
#---------------------------------------------------------------------------------
def getLayers(df, lattice, layerNames):
    
    '''
    Parameters
    ----------
    df
        dataframe : atomic parameters imported from csv file
    lattice
        Lattice : unit cell lattice parameters
    layerNames
        list (str) : layer names as documented in imported csv

    Returns
    -------
    layers : list (Layer)
    '''
    from pyfaults.layerAtom import LayerAtom
    from pyfaults.lattice import Lattice
    
    newLatt = Lattice(a=lattice.a, 
                          b=lattice.b, 
                          c=lattice.c,
                          alpha=lattice.alpha, 
                          beta=lattice.beta, 
                          gamma=lattice.gamma)

    layers = []
    for i in range(len(layerNames)):
        alist = []
        for index, row in df.iterrows():
            if row['Layer'] == layerNames[i]:
                xyz = [row['x'], row['y'], row['z']]
                    
                # create new LayerAtom instance
                newAtom = LayerAtom(layerNames[i], 
                                    row['Atom'], 
                                    row['Element'], 
                                    xyz, 
                                    row['Occupancy'],
                                    row['Biso'],
                                    newLatt)
                alist.append(newAtom)
                
        # create new Layer instance
        newLayer = Layer(alist, newLatt, layerNames[i])
        layers.append(newLayer)
    return layers