'''
Layer class -- single layer of atoms within unit cell lattice
Child_Layer class -- copy of Layer shifted by a translation vector

Layer construction functions:
    - import_CSV
    - get_layer
    - gen_child
'''

import os
import pandas as pd
from diffpy.structure.structure import Atom
import copy as cp
import numpy as np

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' LAYER CLASS '''
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class Layer(list):
    def __init__(self, atoms, lattice, layer_name):
        '''
        Layer class initialization

        Parameters
        ----------
        atoms : list (Atom)
            Atoms in layer, uses Atom object from diffpy
        lattice : Lattice
            Unit cell lattice parameters, uses Lattice object from diffpy
        layer_name : str
            Unique layer name

        Returns
        -------
        None

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
    
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    ''' LAYER CLASS METHODS '''
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    
    def show_labels(self):
        '''
        Prints atom labels and [x,y,z] of all atoms in layer

        Returns
        -------
        None

        '''
        atoms = self.atoms
        for i in atoms:
            print(i.label, i.xyz) 


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' CHILD_LAYER CLASS '''
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class Child_Layer(list):
    def __init__(self, layer_name, parent, shift):
        '''
        Child_Layer class initialization

        Parameters
        ----------
        layer_name : str
            Unique layer name
        parent : Layer
            Parent layer to inheret atoms and lattice from to generate child
        shift : array
            Translation vector [x, y, z]

        Returns
        -------
        None

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
    
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    ''' LAYER CLASS METHODS '''
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    
    def show_labels(self):
        '''
        Prints atom labels and [x,y,z] of all atoms in layer

        Returns
        -------
        None

        '''
        atoms = self.atoms
        for i in atoms:
            print(i.label, i.xyz)
            

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' LAYER CONSTRUCTION FUNCTIONS ''' 
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
    
def import_csv(path, fn):
    '''
    Imports dataframe from csv file

    Parameters
    ----------
    path : str
        File directory
    fn : str
        File name

    Returns
    -------
    df : dataframe
        Data from csv file

    '''

    csv_path = os.path.join(path, fn + ".csv")
    df = pd.read_csv(csv_path)
    return df

#------------------------------------------------------------------------------
def get_layer(df, layer_name, latt):
    '''
    Generates layer from unit cell data
    Assumes each atom entry has parameters Layer, Atom, Element, x, y, z, Occupancy

    Parameters
    ----------
    df : dataframe
        Unit cell information
    layer_name : str
        Unique tag to identify layer
    latt : Lattice
        Unit cell lattice parameters, uses Lattice object from diffpy

    Returns
    -------
    new_layer : Layer
        New layer object

    '''

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
def gen_child(name, parent, t_vect):
    '''
    Generates child layer from a parent

    Parameters
    ----------
    name : str
        Unique layer name
    parent : Layer
        Parent layer
    t_vect : array
        Translation vector [x, y, z]

    Returns
    -------
    child_layer : Child_Layer
        New Child_Layer object

    '''

    child_layer = Child_Layer(layer_name=name, parent=parent, shift=t_vect)
    return child_layer
    
  
    
    
    
    
    
    
    