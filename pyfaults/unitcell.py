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
        for i in range(len(layerList)):
            print(layerList[i].layerName)
    
    # prints atom labels and positions for each layer ----------------------------
    def atom_info(self):
        for i in self.layers:
            for j in i.atoms:
                print(i.layerName + ", " + j.atomLabel + ", " + str(j.xyz))


''' UNIT CELL GENERATOR METHOD ''' 
#---------------------------------------------------------------------------------
# generates new unit cell from CSV -----------------------------------------------
#---------------------------------------------------------------------------------
def genUnitCell(cellName, path, csvFile, lattice, layerList, stackDir, 
                childLayers=None):
    '''
    Parameters
    ----------
    cellName
        str : unique identifier for Unitcell instance
    path
        str : atomic parameters CSV file directory
    csvFile
        str : atomic parameters CSV file name
    lattice
        Lattice : unit cell lattice parameters
    layerList
        list (str) : unique layer names defined in CSV
    stackDir
        str : stacking direction ('a', 'b', or 'c')
    childLayers
        list (str, nparray) [optional] : child layer parameters,
        format each list entry as
        [child name, parent name, translation vector]
        The default is None.

    Returns
    -------
    unitcell : Unitcell
    '''
    from pyfaults import importCSV
    import pyfaults.layer as lyr
    
    # import CSV data
    csv = importCSV(path, csvFile)
    layers = lyr.getLayers(csv, lattice, layerList, stackDir)
    
    # generate child layers if applicable
    if childLayers is not None:
        for i in range(len(childLayers)):
            name = childLayers[i][0]
            parent = childLayers[i][1]
            vec = childLayers[i][2]
            
            for j in layers:
                if j.layerName == parent:
                    childLayer = j.genChildLayer(name, vec)
                    layers.append(childLayer)
                    
    # new Unitcell instance
    unitcell = Unitcell(cellName, layers, lattice)
    return unitcell