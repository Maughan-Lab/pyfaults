##################################################################################
# pyfaults
# Author: Sinclair R. Combs
##################################################################################

import pandas as pd
import numpy as np
import re

#---------------------------------------------------------------------------------
# submodules ---------------------------------------------------------------------
#---------------------------------------------------------------------------------
from pyfaults.lattice import Lattice
from pyfaults.layerAtom import LayerAtom
from pyfaults.layer import Layer
from pyfaults.unitcell import Unitcell
from pyfaults.supercell import Supercell

import pyfaults.layer
import pyfaults.unitcell
import pyfaults.genUnitCell
import pyfaults.supercell
import pyfaults.genSupercells
import pyfaults.simXRD
import pyfaults.calcSims
import pyfaults.plotXRD
import pyfaults.diffCurve
import pyfaults.calcDiffs
import pyfaults.calcFitDiffs
import pyfaults.compareFits
import pyfaults.autoSearch


''' TOP LEVEL METHODS '''
#---------------------------------------------------------------------------------
# writes unit cell or supercell to CIF file --------------------------------------
#---------------------------------------------------------------------------------
def toCif(cell, path, fn):
    '''
    Parameters
    ----------
    cell 
        Unitcell or Supercell : structure to parse to cif file format
    path
        str : file save directory
    fn
        str : file name
    '''
    
    lines = []
    # space group info
    lines.extend([
        '%-31s %s' % ('_symmetry_space_group_name_H-M', 'P1'),
        '%-31s %s' % ('_symmetry_Int_Tables_number', '1'),
        '%-31s %s' % ('_symmetry_cell_setting', 'triclinic'),
        ''])
    
    # lattice parameters
    lines.extend([
        '%-31s %.6g' % ('_cell_length_a', cell.lattice.a),
        '%-31s %.6g' % ('_cell_length_b', cell.lattice.b),
        '%-31s %.6g' % ('_cell_length_c', cell.lattice.c),
        '%-31s %.6g' % ('_cell_angle_alpha', cell.lattice.alpha),
        '%-31s %.6g' % ('_cell_angle_beta', cell.lattice.beta),
        '%-31s %.6g' % ('_cell_angle_gamma', cell.lattice.gamma),
        ''])
    
    # loop info
    lines.extend([
        'loop_',
        '  _atom_site_label',
        '  _atom_site_type_symbol',
        '  _atom_site_fract_x',
        '  _atom_site_fract_y',
        '  _atom_site_fract_z',
        '  _atom_site_U_iso_or_equiv',
        '  _atom_site_adp_type',
        '  _atom_site_occupancy' ])

    # atoms
    for lyr in cell.layers:
        for a in lyr.atoms:
            label = a.atomLabel
            elem = a.element
            x = a.x
            y = a.y
            z = a.z
            occ = a.occupancy
            aline = ' %-5s %-3s %11.6f %11.6f %11.6f %11.6f %-5s %.4f' % (
                label, elem, x, y, z, 2.0, 'Uiso', occ)
            lines.append(aline)

    with open(path + fn + '.cif', 'w') as cif:
        for i in lines:
            cif.write(i + '\n')
    cif.close()

#---------------------------------------------------------------------------------
# converts 2theta values to Q ----------------------------------------------------
#---------------------------------------------------------------------------------
def tt_to_q(twotheta, wavelength):
    '''
    Parameters
    ----------
    twotheta
        nparray : 2theta values (degrees)
    wavelength
        float : instrument wavelength (A)

    Returns
    -------
    Q : nparray
    '''
    Q = 4 * np.pi * np.sin((twotheta * np.pi)/360) / wavelength
    return Q

#---------------------------------------------------------------------------------
# normalizes intensity values ----------------------------------------------------
#---------------------------------------------------------------------------------
def norm(ints):
    '''
    Parameters
    ----------
    ints
        nparray : intensity values

    Returns
    -------
    norm_ints : nparray
    '''
    norm_ints = (ints / np.max(ints))
    return norm_ints

#---------------------------------------------------------------------------------
# imports atomic parameters from CSV file ---------------------------------------- 
#---------------------------------------------------------------------------------       
def importCSV(path, fn):
    '''
    Parameters
    ----------
    path
        str : directory of CSV file
    fn
        str : file name

    Returns
    -------
    df : dataframe
    '''
    df = pd.read_csv(path + fn + '.csv')
    return df

#---------------------------------------------------------------------------------
# imports text file with XRD data ------------------------------------------------
#---------------------------------------------------------------------------------
def importFile(path, fn, ext=None):
    '''
    Parameters
    ----------
    path
        str : directory of text file
    fn
        str : file name
    ext
        str [optional] : file extension if not '.txt'

    Returns
    -------
    q : nparray
    ints : nparray
    '''
    
    if ext is None:
        q, ints = np.loadtxt(path + fn + '.txt', unpack=True, dtype=float)
    if ext is not None:
        q, ints = np.loadtxt(path + fn + ext, unpack=True, dtype=float)
    return q, ints
    
#---------------------------------------------------------------------------------
# import experimental XRD and adjusts 2theta range to match simulated XRD --------
#---------------------------------------------------------------------------------
def importExpt(path, fn, wl, maxTT, ext=None):
    '''
    Parameters
    ----------
    path
        str : experimental data file directory
    fn
        str : file name
    wl
        float : instrument wavelength (A)
    maxTT
        float : maximum 2theta value (degrees)
    ext
        str [optional] : file extension if not '.txt'

    Returns
    -------
    exptQ
        nparray : Q values (A^-1)
    exptNorm : nparray
        nparray : normalized intensity values
    '''
    # import experimental data
    exptTT, exptInts = importFile(path, fn, ext=ext)
    
    # truncate 2theta range according to maxTT
    truncTT = []
    truncInts = []
    for i in range(len(exptTT)):
        if exptTT[i] <= maxTT:
            truncTT.append(exptTT[i])
            truncInts.append(exptInts[i])
    truncTT = np.array(truncTT)
    truncInts = np.array(truncInts)
    
    # convert to Q
    exptQ = tt_to_q(truncTT, wl)
    # normalize intensities
    exptIntsMin = truncInts - np.min(truncInts)
    exptNorm = norm(exptIntsMin)
    
    return exptQ, exptNorm

#---------------------------------------------------------------------------------
# generates stacking fault vector and probability labels -------------------------
#---------------------------------------------------------------------------------
def formatStr(probList, sVecList):
    '''
    Parameters
    ----------
    probList
        list (float) : fault probabilities (0 to 1)
    sVecList
        list (nparray) : stacking vectors ([x, y, z])

    Returns
    -------
    probStrList : list (str)
    sVecStrList : list (str)
    '''
    probStrList = []
    for p in range(len(probList)):
        probStr = re.sub('x', str(int(probList[p]*100)), r'$P = x \%$')
        probStrList.append(probStr)
        
    sVecStrList = []
    for s in range(len(sVecList)):
        txt = r'$\vec{S}_n = \left[ x, y, z \right]$'
        sVecStr = re.sub('n', str(s+1), txt)
        strVars = ['x', 'y', 'z']
        for i in range(2):
            txtSub = re.sub(strVars[i], str(sVecList[s][i]), sVecStr)
            sVecStr = txtSub
        sVecStrList.append(sVecStr)
        
    return probStrList, sVecStrList
 
    
#---------------------------------------------------------------------------------      
# assert imports -----------------------------------------------------------------
#---------------------------------------------------------------------------------
assert Lattice
assert LayerAtom
assert Layer
assert Unitcell
assert Supercell

assert pyfaults.layer
assert pyfaults.unitcell
assert pyfaults.genUnitCell
assert pyfaults.supercell
assert pyfaults.genSupercells
assert pyfaults.simXRD
assert pyfaults.calcSims
assert pyfaults.plotXRD
assert pyfaults.diffCurve
assert pyfaults.calcDiffs
assert pyfaults.calcFitDiffs
assert pyfaults.compareFits
assert pyfaults.autoSearch