#########################################################################################
# pyfaults.unitcell
# Author: Sinclair R. Combs
#########################################################################################

''' Unitcell object class -- Stores unit cell properties '''
#----------------------------------------------------------------------------------------

class Unitcell(object):
    '''
    Parameters
    ----------
    name : str
        Unique identifier for Unitcell object
    layers : list (Layer)
        List of layers of the unit cell
    lattice : Lattice
        Unit cell lattice parameters
    '''
    
    # properties ------------------------------------------------------------------------
    name = property(lambda self: self._name,
                    lambda self, val: self.setParam(name=val))
    
    layers = property(lambda self: self._layers,
                      lambda self, val: self.setParam(layers=val))
    
    lattice = property(lambda self: self._lattice,
                       lambda self, val: self.setParam(lattice=val))
    
    # creates instance of Unitcell object -----------------------------------------------
    def __init__(self, name, layers, lattice):
        
        # initialize parameters
        self._name = None
        self._layers = None
        self._lattice = None
        
        self.setParam(name, layers, lattice)
        
        return
    
    # sets cell parameters --------------------------------------------------------------
    def setParam(self, name=None, layers=None, lattice=None):
        if name is not None:
            self._name = name
        if layers is not None:
            self._layers = layers
        if lattice is not None:
            self._lattice = lattice
            
        return
    
    # prints layer names ----------------------------------------------------------------
    def layer_info(self):
        layerList = self.layers
        for i in range(len(layerList)):
            print(layerList[i].layerName)
    
    # prints atom labels and positions for each layer -----------------------------------
    def atom_info(self):
        for i in self.layers:
            for j in i.atoms:
                print(i.layerName + ", " + j.atomLabel + ", " + str(j.xyz))
