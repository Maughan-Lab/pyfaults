"""
Functions for writing to external files
"""

import os
import copy as cp
import pickle as p

from diffpy.structure.structure import Structure
from diffpy.structure.parsers import p_cif


#------------------------------------------------------------------------------
''' converts a unit cell or supercell object to a structure object
REQUIRED for converting to cif file '''

def write_struct(obj):
    #--------------------------------------------------------------------------
    # obj -- Unitcell, UF_supercell, or FLT_supercell object
    
    # returns Structure object
    #--------------------------------------------------------------------------
    
    # create new structure object
    struct = Structure()
    
    # set structure lattice to be unit cell or supercell lattice
    struct.lattice = obj.lattice
    
    # create a copy of each atom in each layer and add to structure
    layers = obj
    for l in layers:
        for a in l:
            atom = cp.deepcopy(a)
            struct.addNewAtom(atom)
    
    return struct


#------------------------------------------------------------------------------
''' writes a CIF file from a given structure '''

def cif_file(directory, filename, obj):
    #--------------------------------------------------------------------------
    # directory (str) -- file path to save cif to
    # filename (str) -- file name 
    # obj -- Structure object to convert to cif format
    #--------------------------------------------------------------------------
    
    # create .cif file name in given directory
    path = os.path.join(directory, filename + ".cif")
    
    # convert structure object to lines of strings
    obj_list = p_cif.P_cif().toLines(obj)
    
    # create new file
    new_file = open(path, "w")
    
    # write structure information to new cif file
    for i in obj_list:
        new_file.write(i + "\n")
    new_file.close()
    
    
#------------------------------------------------------------------------------
''' write unit cell info to external pickle file '''    
    
def write_pickle(cell, directory, name):
    #--------------------------------------------------------------------------
    # cell -- Unitcell, UF_supercell, or FLT_supercell object
    # directory (str) -- file path to save pickle file to
    # name (str) -- file name
    #--------------------------------------------------------------------------
    
    # create .pickle file path
    path = os.path.join(directory, name + ".pickle")
    
    # create new file in path
    new_file = open(path, "wb")
    
    # write cell information to new file
    p.dump(cell, new_file)
    new_file.close()
    
    
#------------------------------------------------------------------------------
''' read unit cell info from external pickle file '''
    
def load_pickle(directory, name):
    #--------------------------------------------------------------------------
    # directory (str) -- file path to read pickle file from
    # name (str) -- file name
    
    # returns Unitcell, UF_supercell, or FLT_supercell object
    #--------------------------------------------------------------------------
    
    # create .pickle file path
    path = os.path.join(directory, name + ".pickle")
    
    # open file in path
    load_file = open(path, "rb")
    
    # retrieve data from file
    load_data = p.load(load_file)
    load_file.close()
    
    return load_data    
    
    