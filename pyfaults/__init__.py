#########################################################################################
# pyfaults
# Author: Sinclair R. Combs
#########################################################################################

''' 
Classes:
    Lattice
    LayerAtom
    Layer
    Child
    Unitcell
    Supercell
'''

import pandas as pd

# submodules ----------------------------------------------------------------------------

from pyfaults.lattice import Lattice
from pyfaults.layerAtom import LayerAtom
from pyfaults.layer import Layer
from pyfaults.childLayer import ChildLayer
from pyfaults.unitcell import Unitcell
from pyfaults.supercell import Supercell

import pyfaults.simXRD

# top level methods ---------------------------------------------------------------------

def toCif(cell, path, fn):
    '''
    Writes cif file from unit cell or supercell object

    Parameters
    ----------
    cell : Unitcell or Supercell
        Structure to parse to cif file format
    path : str
        File directory to save cif to
    fn : str
        File name
    '''
    
    lines = []
    
    # space group info
    lines.extend([
        "%-31s %s" % ("_symmetry_space_group_name_H-M", "'P1"),
        "%-31s %s" % ("_symmetry_Int_Tables_number", "1"),
        "%-31s %s" % ("_symmetry_cell_setting", "triclinic"),
        ""])
    
    # lattice parameters
    lines.extend([
        "%-31s %.6g" % ("_cell_length_a", cell.lattice.a),
        "%-31s %.6g" % ("_cell_length_b", cell.lattice.b),
        "%-31s %.6g" % ("_cell_length_c", cell.lattice.c),
        "%-31s %.6g" % ("_cell_angle_alpha", cell.lattice.alpha),
        "%-31s %.6g" % ("_cell_angle_beta", cell.lattice.beta),
        "%-31s %.6g" % ("_cell_angle_gamma", cell.lattice.gamma),
        ""])
    
    # loop info
    lines.extend([
        "loop_",
        "  _atom_site_label",
        "  _atom_site_type_symbol",
        "  _atom_site_fract_x",
        "  _atom_site_fract_y",
        "  _atom_site_fract_z",
        "  _atom_site_U_iso_or_equiv",
        "  _atom_site_adp_type",
        "  _atom_site_occupancy" ])

    # atoms
    for lyr in cell.layers:
        for a in lyr.atoms:
            label = a.atomLabel
            elem = a.element
            x = a.x
            y = a.y
            z = a.z
            occ = a.occupancy
            aline = " %-5s %-3s %11.6f %11.6f %11.6f %11.6f %-5s %.4f" % (
                label, elem, x, y, z, 2.0, "Uiso", occ)
            lines.append(aline)

    with open(path + fn + ".cif", "w") as cif:
        for i in lines:
            cif.write(i + "\n")
    cif.close()
        
        
def importCSV(path, fn):
    '''
    Imports atomic parameters from csv file
    
    Parameters
    ----------
    path : str
        File directory where csv is stored
    fn : str
        File name

    Returns
    -------
    df : dataframe
        Atomic parameters
    '''
    df = pd.read_csv(path + fn + ".csv")
    return df
        
        
assert Lattice
assert LayerAtom
assert Layer
assert ChildLayer
assert Unitcell
assert Supercell
assert pyfaults.simXRD