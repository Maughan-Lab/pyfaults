#########################################################################################
# pyfaults.layerAtom
# Author: Sinclair R. Combs
#########################################################################################

''' LayerAtom object class -- Stores properties of atom in specific layer '''
#----------------------------------------------------------------------------------------

class LayerAtom(object):
    '''
    Parameters
    ----------
    layerName : str
        Unique identifier for layer in which atom resides
    atomLabel : str
        Unique identifier for LayerAtom object
    element : str
        Atomic element and oxidation state if applicable
    xyz : nparray
        Atomic position in fractional coordinates
    occupancy : float
        Crystallographic site occupancy
    lattice : Lattice
        Unit cell lattice parameters
    '''
        
    # properties ------------------------------------------------------------------------
    layerName = property(lambda self: self._layerName,
                         lambda self, val: self.setParam(layerName=val))
        
    atomLabel = property(lambda self: self._atomLabel,
                         lambda self, val: self.setParam(atomLabel=val))
        
    element = property(lambda self: self._element,
                       lambda self, val: self.setParam(element=val))
        
    xyz = property(lambda self: self._xyz,
                   lambda self, val: self.setParam(xyz=val))
        
    x = property(lambda self: self._xyz[0],
                 lambda self, val: self.setParam(0, val),
                 doc="float : fractional coordinate x")
        
    y = property(lambda self: self._xyz[1],
                 lambda self, val: self.setParam(1, val),
                 doc="float : fractional coordinate y")
        
    z = property(lambda self: self._xyz[2],
                 lambda self, val: self.setParam(2, val),
                 doc="float : fractional coordinate z")
        
    occupancy = property(lambda self: self._occupancy,
                         lambda self, val: self.setParam(occupancy=val))
        
    lattice = property(lambda self: self._lattice,
                       lambda self, val: self.setParam(lattice=val))
    
    # creates instance of LayerAtom object ----------------------------------------------
    def __init__(self, layerName, atomLabel, element, xyz, occupancy, lattice):
        
        # initialize parameters
        self._layerName = None
        self._atomLabel = None
        self._element = None
        self._xyz = None
        self._x = None
        self._y = None
        self._z = None
        self._lattice = None
        self._occupancy = None
        
        self.setParam(layerName, atomLabel, element, xyz, lattice, occupancy)
        
        return
    
    # sets atom parameters --------------------------------------------------------------
    def setParam(self, layerName=None, atomLabel=None, element=None, xyz=None,
                 lattice=None, occupancy=None):
        
        if layerName is not None:
            self._layerName = layerName
        if atomLabel is not None:
            self._atomLabel = atomLabel + "_" + layerName
        if element is not None:
            self._element = element
        if xyz is not None:
            self._xyz = xyz
            self._x = xyz[0]
            self._y = xyz[1]
            self._z = xyz[2]
        if lattice is not None:
            self._lattice = lattice
        if occupancy is not None:
            self._occupancy = occupancy
        
        return
    
    # prints atom layer, name, element, position, and occupancy -------------------------
    def display(self):
        a = "\n".join(("Layer Name: " + self.layerName, 
                       "Atom Label: " + self.atomLabel,
                       "Element: " + self.element,
                       "Atomic Position: " + str(self.xyz),
                       "Occupancy: " + str(self.occupancy)))
        print(a)