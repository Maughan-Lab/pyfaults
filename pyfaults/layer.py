##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import copy as cp
import numpy as np

# Layer object class ----------
class Layer(object):
'''
Parameters
----------
atoms (array_like of LayerAtom) : list of atoms in layer
lattice (Lattice) : unit cell lattice parameters
layerName (str) : unique identifier for layer
'''   
    # properties ----------
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
    
    # initialization, defines Layer defaults ----------
    def __init__(self, atoms, lattice, layerName):
        self._layerName = None
        self._atoms = None
        self._lattice = None
        self.setParam(atoms, lattice, layerName)
        return
    
    # sets parameters of Layer ----------
    def setParam(self, atoms=None, lattice=None, layerName=None):
        if atoms is not None:
            self._atoms = atoms
        if lattice is not None:
            self._lattice = lattice
        if layerName is not None:
            self._layerName = layerName
        return
    
    # prints atoms in Layer ----------
    def display(self):
        print(self.layerName)
        for i in self.atoms:
            print(i.atomLabel + ', ' + str(i.xyz))
        return 
            
    # generates CIF file of layer ----------
    def toCif(self, path):
    '''
    Parameters
    ----------
    path (str) : file path
    '''
        from pyfaults import toCif
        toCif(self, path, 'Layer_' + self.layerName + '_CIF')
    
    # creates a new Layer object that is a copy of the original Layer displaced by a given translation vector ----------
    def genChildLayer(self, childName, transVec):
        '''
        Parameters
        ----------
        childName (str) : unique identifier for layer
        transVec (array_like) : translation vector [x, y, z] relative to original/parent layer

        Returns
        -------
        childLayer (Layer) : new child layer
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

# imports layer information from a DataFrame ----------
def getLayers(df, lattice, layerNames):
    
    '''
    Parameters
    ----------
    df (DataFrame) : atomic parameters
    lattice (Lattice) : unit cell lattice parameters
    layerNames (array_like) : defined layer names

    Returns
    -------
    layers (array_like) : list of Layer objects
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
