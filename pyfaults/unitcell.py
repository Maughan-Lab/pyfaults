'''
Defines Unitcell class
'''

import copy as cp


#------------------------------------------------------------------------------
''' Unitcell class initialization '''

class Unitcell(list):
    def __init__(self, layers):
        #----------------------------------------------------------------------
        # layers -- list of Layer instances
        #----------------------------------------------------------------------
        
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

#------------------------------------------------------------------------------ 
# Unitcell methods
    
    #--------------------------------------------------------------------------
    ''' prints atom labels and [x,y,z] of all atoms in unit cell'''
    
    def show_atoms(self):
        for l in self:
            for a in l:
                print(a.label, a.xyz)
    
    
    #--------------------------------------------------------------------------
    '''' prints name of all layers in unit cell'''
    
    def show_layers(self):
        layer = self.layers
        for i in layer:
            i.layer_name
        
        
        
        