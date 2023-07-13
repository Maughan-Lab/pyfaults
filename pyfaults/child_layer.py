'''
Defines Child_Layer class
'''

import copy as cp
import numpy as np

class Child_Layer(list):
    def __init__(self, layer_name, parent, shift):
        '''
        layer_name (str)
        parent (Layer)
        shift (nparray)
        '''
        
        # create and assign lattice variable based on parent layer lattice
        self.lattice = parent.lattice
        
        # assign layer_name, parent, and shift variables
        self.layer_name = layer_name
        self.parent = parent
        self.shift = shift
        
        # create a copy of parent layer atoms
        child_atoms = cp.deepcopy(parent.atoms)
        
        # add translation vector and layer name to each child atom
        for a in child_atoms:
            a.xyz = np.add(a.xyz, shift)
            label = a.label.split("_")
            a.label = label[0] + "_" + layer_name
        
        # assign atoms variable
        self.atoms = child_atoms
        
        # appends atoms to the Child_Layer object (list)
        for a in child_atoms:
            self.append(a)

        return

#------------------------------------------------------------------------------ 
    ''' Child_Layer methods '''

    # prints atom labels and [x,y,z] coordinates
    def show_labels(self):
        atoms = self.atoms
        for i in atoms:
            print(i.label, i.xyz)
