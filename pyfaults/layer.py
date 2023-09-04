#########################################################################################
# pyfaults.layer
# Author: Sinclair R. Combs
#########################################################################################

import copy as cp
import numpy as np

''' Layer object class -- Stores layer properties'''
#----------------------------------------------------------------------------------------

class Layer(object):
    '''
    Parameters
    ----------
    atoms : list (LayerAtom)
        Atoms in layer
    lattice : Lattice
        Unit cell lattice parameters
    layerName : str
        Unique identifier for Layer object
    '''
    
    # properties ------------------------------------------------------------------------
    atoms = property(lambda self: self._atoms,
                     lambda self, val: self.setParam(atoms=val))
    
    lattice = property(lambda self: self._lattice,
                       lambda self, val: self.setParam(lattice=val))
    
    layerName = property(lambda self: self._layerName,
                         lambda self, val: self.setParam(layerName=val))
    
    # creates instance of Layer object --------------------------------------------------
    def __init__(self, atoms, lattice, layerName):
        
        # initialize parameters
        self._layerName = None
        self._atoms = None
        self._lattice = None

        self.setParam(atoms, lattice, layerName)

        return
    
    # sets layer parameters -------------------------------------------------------------
    def setParam(self, atoms=None, lattice=None, layerName=None):
        if atoms is not None:
            self._atoms = atoms
        if lattice is not None:
            self._lattice = lattice
        if layerName is not None:
            self._layerName = layerName
        
        return
    
    # prints layer name and atoms -------------------------------------------------------
    def display(self):
        lyr = "\n".join(("Layer Name: " + self.layerName, 
                       "Atoms: " + str(self.atoms)))
        print(lyr)
    
    # prints atom labels and positions ---------------------------------------------------
    def atom_info(self):
        for i in self.atoms:
            print(i.atomLabel + ", " + str(i.xyz))


''' Pulls layers from imported dataframe'''
#----------------------------------------------------------------------------------------
            
def getLayers(df, lattice, layerNames, stackDir):
    from pyfaults.layerAtom import LayerAtom
    from pyfaults.lattice import Lattice
    
    '''
    Parameters
    ----------
    df : dataframe
        Atomic parameters imported from csv file
    lattice : Lattice
        Unit cell lattice parameters
    layerNames : list (str)
        Layer names as documented in imported csv
    stackDir : str
        Lattice vector orthogonal to layers

    Returns
    -------
    layers : list (Layer)
    '''
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
    
    return layers

''' Generates a child layer '''
#----------------------------------------------------------------------------------------
def genChildLayer(self, childName, transVec):
    from pyfaults.layerAtom import LayerAtom
    
    childAtoms = []
    
    for a in self.atoms:
        pAtom = cp.deepcopy(a)
        
        splitLabel = pAtom.atomLabel.split("_")
        newLabel = splitLabel[0] + "_" + childName
        
        newPos = np.add(pAtom.xyz, transVec)
        
        cAtom = LayerAtom(childName, newLabel, pAtom.element, newPos,
                          pAtom.occupancy, self.lattice)
        
        childAtoms.append(cAtom)
        
    childLayer = Layer(childAtoms, self.lattice, childName)
    
    return childLayer






