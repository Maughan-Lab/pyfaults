'''
Defines UF_supercell class
'''

import copy as cp

class UF_supercell(list):
    def __init__(self, unitcell, n_stacks):
        '''
        unitcell (UF_unitcell)
        n_stacks (int)
        '''
        
        # assign n_stacks variable
        self.n_stacks = n_stacks
        
        # create a copy of the unit cell lattice
        latt = cp.deepcopy(unitcell.lattice)
        
        # scale the c-axis vector by the number of stacks
        latt.setLatPar(c = n_stacks * latt.c)
        
        # create and assign lattice variable
        self.lattice = latt
        
        # create a copy of the unit cell for each N in n_stacks
        for n in range(0, n_stacks):
            uc = cp.deepcopy(unitcell)
            
            # create a copy of each layer in the unit cell
            for layer in uc:
                wl = cp.deepcopy(layer)
                
                # scale the z-position of each atom and add label based on N
                for a in wl:
                    a.z = (a.z + n) / n_stacks
                    a.label = a.label + "_n" + str(n)
                
                # set lattice of layer to N-scaled lattice
                wl.lattice = latt
                
                # add label to layer name based on N
                wl.layer_name = wl.layer_name + "_n" + str(n)
                
                # append layers to the UF_supercell object (list)
                self.append(wl)

        return

#------------------------------------------------------------------------------ 
    ''' UF_supercell methods '''
    
    # prints atom labels and [x,y,z] coordinates
    def show_atoms(self):
        for l in self:
            for a in l:
                print(a.label, a.xyz)
    
    # prints layer names
    def show_layers(self):
        layer = self.layers
        for i in layer:
            i.layer_name

