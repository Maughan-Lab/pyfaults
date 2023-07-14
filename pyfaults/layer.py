'''
Defines Layer and Child_Layer classes
'''

import os
import pandas as pd
from diffpy.structure.structure import Atom
import copy as cp
import numpy as np


#------------------------------------------------------------------------------
''' Layer class initialization '''

class Layer(list):
    def __init__(self, atoms, lattice, layer_name):
        #----------------------------------------------------------------------
        # atoms -- list of Atom objects
        # lattice -- Lattice object
        # layer_name (str) -- unique layer name
        #----------------------------------------------------------------------
        
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
# Layer methods
    
    #--------------------------------------------------------------------------
    ''' prints atom labels and [x,y,z] of all atoms in layer'''
    
    def show_labels(self):
        atoms = self.atoms
        for i in atoms:
            print(i.label, i.xyz) 


#------------------------------------------------------------------------------
''' Child_Layer class initialization '''

class Child_Layer(list):
    def __init__(self, layer_name, parent, shift):
        #----------------------------------------------------------------------
        # layer_name (str) -- unique layer name
        # parent -- parent Layer object
        # shift (nparray) -- translation vector [x, y, z]
        #----------------------------------------------------------------------
        
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
# Child_Layer methods

    #--------------------------------------------------------------------------
    ''' prints atom labels and [x,y,z] of all atoms in layer'''
    
    def show_labels(self):
        atoms = self.atoms
        for i in atoms:
            print(i.label, i.xyz)
            

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' Functions for generating layers from unit cell csv files ''' 
    

#------------------------------------------------------------------------------
''' import dataframe from csv file '''

def import_csv(directory, filename):
    #--------------------------------------------------------------------------
    # directory (str) -- file path to read csv file from
    # name (str) -- file name
    
    # returns dataframe of unit cell information
    #--------------------------------------------------------------------------
    
    # create .csv file path
    path = os.path.join(directory, filename + ".csv")
    
    # retrieve data from file
    df = pd.read_csv(path)
    
    return df


#------------------------------------------------------------------------------
''' generate layer from csv data '''
    
def get_layer(df, layer_name, latt):
    #--------------------------------------------------------------------------
    # df -- dataframe containing unit cell information
    # layer_name (str) -- unique layer tag in csv data
    # latt -- Lattice object for unit cell
    
    # returns Layer object
    #--------------------------------------------------------------------------
    
    # retrieve atoms corresponding to given layer name
    l = df[df["Layer"] == layer_name]
    
    # create new Atom objects for each atom in layer
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
    
    # create new Layer object
    new_layer = Layer(atom_list, latt, layer_name)
    
    return new_layer


#------------------------------------------------------------------------------
''' generate child layer from parent layer (basically creates new Child_layer
instance) '''

def gen_child(name, parent, t_vect):
    #--------------------------------------------------------------------------
    # name (str) -- name of child layer
    #parent = parent Layer object
    # t_vect (nparray) -- translation vector [x, y, z]
    
    # returns Child_Layer object
    #--------------------------------------------------------------------------
    child_layer = Child_Layer(layer_name=name, parent=parent, shift=t_vect)
    
    return child_layer
    
  
    
    
    
    
    
    
    