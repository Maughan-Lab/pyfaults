#########################################################################################
# Author: Sinclair R. Combs
#########################################################################################

''' UNITCELL CLASS '''
#---------------------------------------------------------------------------------
class Unitcell(object):
    # properties -----------------------------------------------------------------
    name =\
        property(lambda self: self._name,
                 lambda self, val: self.setParam(name=val),
                 doc='str : unique identifier for Unitcell instance')
    
    layers =\
        property(lambda self: self._layers,
                 lambda self, val: self.setParam(layers=val),
                 doc='list (Layer) : layers contained in unit cell')
    
    lattice =\
        property(lambda self: self._lattice,
                 lambda self, val: self.setParam(lattice=val),
                 doc='Lattice : unit cell lattice parameters')
    
    # creates instance of Unitcell object ----------------------------------------
    def __init__(self, name, layers, lattice):
        # initialize parameters
        self._name = None
        self._layers = None
        self._lattice = None
        self.setParam(name, layers, lattice)
        return
    
    # sets parameters -------------------------------------------------------------
    def setParam(self, name=None, layers=None, lattice=None):
        if name is not None:
            self._name = name
        if layers is not None:
            self._layers = layers
        if lattice is not None:
            self._lattice = lattice
        return
    
    # prints layer names ---------------------------------------------------------
    def layer_info(self):
        layerList = self.layers
        layerStr = []
        for i in range(len(layerList)):
            layerStr.append(layerList[i].layerName)
            print(layerList[i].layerName)
        return
                
    # generates CIF of layer -----------------------------------------------------
    def toCif(self, path):
        from pyfaults import toCif
        toCif(self, path, 'Unitcell_' + self.name + '_CIF')