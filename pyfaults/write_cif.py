"""
Functions for writing to external files
"""

import os
import copy as cp

from diffpy.structure.structure import Structure
from diffpy.structure.parsers import p_cif

#------------------------------------------------------------------------------

'''
converts a unit cell or supercell object to a structure object (needed to write 
to a CIF file)
'''
def write_struct(obj):
    '''
    obj (UF_unitcell, FLT_unitcell, UF_supercell, or FLT_supercell)
    '''
    
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

''' writes a CIF file given a structure '''
def cif_file(directory, filename, obj):
    '''
    directory (str)
    filename (str)
    obj (Structure)
    '''
    
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