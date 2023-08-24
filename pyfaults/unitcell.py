'''
Unitcell class -- collection of Layer objects constructing a unit cell
'''

import copy as cp

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' UNITCELL CLASS '''
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class Unitcell(list):
    def __init__(self, layers):
        '''
        Unitcell class initialization

        Parameters
        ----------
        layers : list (Layer)
            List of layer to include in unit cell

        Returns
        -------
        None

        '''

        # create and assign lattice variable based on parent layer lattice
        self.lattice = layers[0].lattice
        
        # create a copy of layers
        uc_layers = cp.deepcopy(layers)
        
        # assign layers variable
        self.layers = uc_layers
        
        # appends layers to the UF_unitcell object (list)
        for l in uc_layers:
            self.append(l)
        
        return
    
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    ''' UNITCELL CLASS METHODS '''
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------

    def show_atoms(self):
        '''
        Prints atom labels and [x,y,z] of all atoms in layer

        Returns
        -------
        None

        '''
        for l in self:
            for a in l:
                print(a.label, a.xyz)
    
    #--------------------------------------------------------------------------
    def show_layers(self):
        '''
        Prints labels all layers in unit cell

        Returns
        -------
        None

        '''
        layer = self.layers
        for i in layer:
            i.layer_name
        
        
        
        