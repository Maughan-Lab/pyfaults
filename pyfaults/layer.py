'''
Defines Layer and Child_Layer classes
'''

import os
import pandas as pd
from diffpy.structure.structure import Atom
import copy as cp
import numpy as np


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
#------------------------------------------------------------------------------ 
    ''' Layer methods '''
    
    # prints atom labels and [x,y,z] coordinates
    def show_labels(self):
        atoms = self.atoms
        for i in atoms:
            print(i.label, i.xyz) 


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
    
#------------------------------------------------------------------------------ 
    ''' Child_Layer methods '''

    # prints atom labels and [x,y,z] coordinates
    def show_labels(self):
        atoms = self.atoms
        for i in atoms:
            print(i.label, i.xyz)
            

#------------------------------------------------------------------------------ 
    
# import dataframe from csv file
def import_csv(directory, filename):
    path = os.path.join(directory, filename + ".csv")
    
    df = pd.read_csv(path)
    
    return df
    
# generate layer from csv data
def get_layer(df, layer_name, latt):
    
    l = df[df["Layer"] == layer_name]
    
    atom_list = []
    for j in range(0, len(l.index)):
        atom = l.iloc[j]["Atom"]
        elem = l.iloc[j]["Element"]
        x = l.iloc[j]["x"]
        y = l.iloc[j]["y"]
        z = l.iloc[j]["z"]
        occ = l.iloc[j]["Occupancy"]
        
        pos = [x,y,z]
        
        new_atom = Atom(atype=elem, xyz=pos, label=atom, occupancy=occ, lattice=latt)
        
        atom_list.append(new_atom)
        
    new_layer = Layer(atom_list, latt, layer_name)
    
    return new_layer
    
# generate child layer
def gen_child(name, parent_layer, t_vect):
    child_layer = Child_Layer(layer_name=name, parent=parent_layer, shift=t_vect)
    
    return child_layer
    
  
    
    
    
    
    
    
    