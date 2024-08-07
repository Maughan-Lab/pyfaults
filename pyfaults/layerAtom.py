##################################################################################
# Author: Sinclair R. Combs
##################################################################################

# LayerAtom object class ----------
class LayerAtom(object):
    '''
    Parameters
    ----------
    layerName (str) -- name of layer that atom resides in
    atomLabel (str) -- unique identifier for atomic position
    element (str) -- elemental species and oxidation state
    xyz (array_like) -- atomic position in fractional coordinates
    occupancy (float) -- site occupancy (0 to 1)
    biso (float) -- isotropic atomic displacement parameter
    lattice (Lattice) -- unit cell lattice parameters
    '''
    # properties ----------
    layerName =\
        property(lambda self: self._layerName,
                 lambda self, val: self.setParam(layerName=val),
                 doc='str : unique identifier for resident layer')
        
    atomLabel =\
        property(lambda self: self._atomLabel,
                 lambda self, val: self.setParam(atomLabel=val),
                 doc='str : unique identifier for LayerAtom instance')
        
    element =\
        property(lambda self: self._element,
                 lambda self, val: self.setParam(element=val),
                 doc='str : atomic element and oxidation state')
        
    xyz =\
        property(lambda self: self._xyz,
                 lambda self, val: self.setParam(xyz=val),
                 doc='nparray : atomic position in fractional coordinates')
        
    x =\
        property(lambda self: self._xyz[0],
                 lambda self, val: self.setParam(0, val),
                 doc='float : x-component of atomic position')
        
    y =\
        property(lambda self: self._xyz[1],
                 lambda self, val: self.setParam(1, val),
                 doc='float : y-component of atomic position')
        
    z =\
        property(lambda self: self._xyz[2],
                 lambda self, val: self.setParam(2, val),
                 doc='float : z-component of atomic position')
        
    occupancy =\
        property(lambda self: self._occupancy,
                 lambda self, val: self.setParam(occupancy=val),
                 doc='float : site occupancy, values from 0 to 1')
        
    biso =\
        property(lambda self: self._biso,
                 lambda self, val: self.setParam(biso=val),
                 doc='float : atomic parameter displacement B_iso')
        
    lattice =\
        property(lambda self: self._lattice,
                 lambda self, val: self.setParam(lattice=val),
                 doc='Lattice : unit cell lattice parameters')
    
    # initialization, defines LayerAtom defaults ----------
    def __init__(self, layerName, atomLabel, element, xyz, occupancy, biso, lattice):
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
        self._biso = None
        self.setParam(layerName, atomLabel, element, xyz, lattice, occupancy, biso)
        return
    
    # sets parameters of LayerAtom ----------
    def setParam(self, layerName=None, atomLabel=None, element=None, xyz=None, 
                 lattice=None, occupancy=None, biso=None):
        if layerName is not None:
            self._layerName = layerName
        if atomLabel is not None:
            self._atomLabel = atomLabel + '_' + layerName
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
        if biso is not None:
            self._biso = biso
        
        return
