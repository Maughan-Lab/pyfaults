'''
Functions to write cif and pickle files
'''

import copy as cp
import pickle as p

from diffpy.structure.structure import Structure
from diffpy.structure.parsers import p_cif

#------------------------------------------------------------------------------
def gen_cif(cell, path, fn, inclPkl=False):
    '''
    Generates cif file

    Parameters
    ----------
    cell : Unitcell, UF, or FLT
        Structure to export, can be Unitcell, UF (supercell), or FLT (supercell) object
    path : str
        File directory
    fn : str
        File name
    inclPkl : bool, optional
        If True, generates pickle file with cell object information

    Returns
    -------
    None

    '''
    
    # create new structure object
    struct = Structure()
    struct.lattice = cell.lattice
    
    # copy each atom in layer and add to structure
    layers = cell
    for l in layers:
        for a in l:
            atom = cp.deepcopy(a)
            struct.addNewAtom(atom)
            
    # structure to string
    cell_str = p_cif.P_cif().toLines(cell)
    
    # write cif file
    cif_file = open(path + fn + ".cif", "wb")
    for i in cell_str:
        cif_file.write(i + "\n")
    cif_file.close()
    
    # write pickle file
    if inclPkl == True:
        pkl_file = open(path + fn + ".pkl", "wb")
        p.dump(cell, pkl_file)
        pkl_file.close()
    
#------------------------------------------------------------------------------
def load_pickle(path, fn):
    '''
    Read pickle file

    Parameters
    ----------
    path : str
        File directory
    fn : str
        File name

    Returns
    -------
    data : Unitcell, UF, or FLT
        Retrieved object data

    '''
    file = open(path + fn + ".pkl", "rb")
    data = p.load(file)
    file.close()
    return data