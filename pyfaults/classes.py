"""
Defines classes Layer, Child_Layer, Unitcell, UF_supercell, FLT_supercell
"""

import numpy as np
import copy as cp
import random as r


''' LAYER CLASS '''
#------------------------------------------------------------------------------
class Layer(list):
    def __init__(self, atoms, lattice, layer_name):
        '''
        atoms (list of Atom instances)
        lattice (Lattice)
        layer_name(str)
        '''
        
        # assign lattice and layer_name variables
        self.lattice = lattice
        self.layer_name = layer_name
        
        # add layer name to each atom
        for a in atoms:
            a.label = a.label + "_" + layer_name
        
        # assign atoms variable
        self.atoms = atoms
        
        # appends atoms to the Layer object (list)
        for a in atoms:
            self.append(a)

        return
    
    def show_labels(self):
        atoms = self.atoms
        for i in atoms:
            print(i.label, i.xyz)
            
            
            
''' CHILD LAYER CLASS '''
#------------------------------------------------------------------------------
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

    def show_labels(self):
        atoms = self.atoms
        for i in atoms:
            print(i.label, i.xyz)



''' UNIT CELL CLASS '''
#------------------------------------------------------------------------------
class Unitcell(list):
    def __init__(self, layers):
        ''' 
        layers (list of Layer instances)
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

    def show_atoms(self):
        for l in self:
            for a in l:
                print(a.label, a.xyz)
                
    def show_layers(self):
        layer = self.layers
        for i in layer:
            i.layer_name


''' IDEAL SUPERCELL CLASS'''
#------------------------------------------------------------------------------
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

    def show_atoms(self):
        for l in self:
            for a in l:
                print(a.label, a.xyz)
                
    def show_layers(self):
        layer = self.layers
        for i in layer:
            i.layer_name
            

''' FAULTED SUPERCELL CLASS '''
#------------------------------------------------------------------------------
class FLT_supercell(list):
    def __init__(self, unitcell, n_stacks, flt_layer, shift, p):
        '''
        unitcell (UF_unitcell)
        n_stacks (int)
        flt_layer (str)
        shift (nparray)
        p (float)
        '''
        
        # assign n_stacks, shift, and p variables
        self.n_stacks = n_stacks
        self.shift = shift
        self.p = p
        
        # create a copy of the unit cell lattice
        latt = cp.deepcopy(unitcell.lattice)
        
        # scale the c-axis vector by the number of stacks
        latt.setLatPar(c = n_stacks * latt.c)
        
        # create and assign lattice variable
        self.lattice = latt
        
        # convert decimal probability to a percentage
        percent = p * 100
        
        # create a copy of the unit cell for each N in n_stacks
        for n in range(0, n_stacks):
            n_label = "_n" + str(n)
            flt_cell = cp.deepcopy(unitcell)
            
            for layer in flt_cell:
                # generates a random integer to assess probability of fault
                gen_prob = r.randint(0, 100)
                
                # if layer matches the name of the given faulted layer and the
                # probability condition is met
                if layer.layer_name == flt_layer and gen_prob <= percent:
                    
                    # create a copy of faulted layer
                    flt = cp.deepcopy(layer)
                    
                    # add stacking vector, scale z-position, and add labels to
                    # each atom in the faulted layer
                    for a in flt:
                        a.xyz = np.add(a.xyz, shift)
                        a.z = (a.z + n) / n_stacks
                        a.label = a.label + "_flt" + n_label
                        
                    # set lattice of layer to N-scaled lattice
                    flt.lattice = latt
                    
                    # add label to layer name based on N
                    flt.layer_name = flt.layer_name + "_flt" + n_label
                    
                    # append faulted layer to the FLT_supercell object (list)
                    self.append(flt)
                    
                # if layer matches the name of the given faulted layer and the
                # probability condition is NOT met
                elif layer.layer_name == flt_layer and gen_prob > percent:
                    
                    # create a copy of layer
                    flt = cp.deepcopy(layer)
                    
                    # scale z-position and add labels to each atom in layer
                    for a in flt:
                        a.z = (a.z + n) / n_stacks
                        a.label = a.label + n_label
                        
                    # set lattice of layer to N-scaled lattice
                    flt.lattice = latt
                    
                    # add label to layer name based on N
                    flt.layer_name = flt.layer_name + n_label
                    
                    # append layer to the FLT_supercell object (list)
                    self.append(flt)
                    
                # if layer does not match the name of the given faulted layer
                elif layer.layer_name != flt_layer:
                    
                    # create a copy of layer
                    l = cp.deepcopy(layer)
                    
                    # scale z-position and add labels to each atom in layer
                    for a in l:
                        a.z = (a.z + n) / n_stacks
                        a.label = a.label + n_label
                        
                    # set lattice of layer to N-scaled lattice
                    l.lattice = latt
                    
                    # add label to layer name based on N
                    l.layer_name = l.layer_name + n_label
                    
                    # append layer to the FLT_supercell object (list)
                    self.append(l)

        return
    
    def show_atoms(self):
        for l in self:
            for a in l:
                print(a.label, a.xyz)
                
    def show_layers(self):
        layer = self.layers
        for i in layer:
            i.layer_name