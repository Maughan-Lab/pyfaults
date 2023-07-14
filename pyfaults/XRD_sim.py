"""
Functions for XRD simulations
"""

import os
import Dans_Diffraction as df
import numpy as np

#------------------------------------------------------------------------------
''' simulates X-ray diffraction patterns '''

def full_XRD_sim(cif_dir, cif, wavelength, tt_max, pw, bg):
    #--------------------------------------------------------------------------
    # cif_dir (str) - file path to load cif from
    # cif (str) -- name of cif file
    # wavelength (float) -- instrument wavelength
    # tt_max (float) -- maximum 2theta value
    # pw (float) -- peak width, units of A^-1
    # bg (float) -- average of normal background
    
    # returns lists of calcualted Q and intensity values
    #--------------------------------------------------------------------------
    
    # load cif file as structure 
    path = os.path.join(cif_dir, cif + ".cif")
    struct = df.Crystal(path)
    
    # setup diffraction parameters
    energy_kev = df.fc.wave2energy(wavelength)
    struct.Scatter.setup_scatter("xray")
    wavevector_max = df.fc.calqmag(tt_max, energy_kev)
    
    # simultate XRD pattern
    q, ints = struct.Scatter.generate_powder(wavevector_max, 
                                             peak_width=pw, background=bg,
                                             powder_average=True)
    return q, ints


#-----------------------------------------------------------------------------
''' saves text file of simulated XRD pattern '''

def save_sim(directory, name, q, ints):
    #--------------------------------------------------------------------------
    # directory (str) -- file path to save data to
    # name (str) -- file name 
    # q (list) -- list of calculated Q values
    # ints (list) -- list of calculated intensity values
    #--------------------------------------------------------------------------
    
    # create file path
    path = os.path.join(directory, name + ".txt")
    
    # create new file
    new_file = open(path, "w")
    
    # write simulation data to new text file
    with new_file as f:
        for (q, ints) in zip(q, ints):
            f.write("{0} {1}\n".format(q, ints))
    new_file.close()
    
    
#-----------------------------------------------------------------------------
''' import XRD simulation text file '''

def import_sim(directory, name):
    #--------------------------------------------------------------------------
    # path (str) -- path to folder where text file is stored
    # name (str) -- text file name
    
    # returns lists of Q and intensity values from text file
    #--------------------------------------------------------------------------
    
    # create data file path 
    path = os.path.join(directory, name + ".txt")
    
    # load data from text file
    q, ints = np.loadtxt(path, unpack=True, dtype=float)
    
    return q, ints
    

#-----------------------------------------------------------------------------
''' load cif as structure '''

def load_cif(directory, name):
    #--------------------------------------------------------------------------
    # directory (str) -- path to folder where cif is stored
    # name (str) -- cif file name
    
    # returns structure as Crystal object
    #--------------------------------------------------------------------------
    
    # creates cif file path
    path = os.path.join(directory, name + ".cif")
    
    # convert cif to Crystal object
    struct = df.Crystal(path)
    
    return struct


#-----------------------------------------------------------------------------
''' set up diffraction simulation conditions '''

def sim_setup(struct, wl, max_tt):
    #--------------------------------------------------------------------------
    # struct -- Crystal structure object
    # wl (float) -- instrument wavelength
    # max_tt (float) -- maximum 2theta value
    
    # returns maximum wavevector for simulation
    #--------------------------------------------------------------------------
    
    energy_kev = df.fc.wave2energy(wl)
    struct.Scatter.setup_scatter("xray")
    max_wavevector = df.fc.calqmag(max_tt, energy_kev)
    
    return max_wavevector
    

#-----------------------------------------------------------------------------
''' simulate powder diffraction '''

def sim(struct, max_wv, pw, bg):
    #--------------------------------------------------------------------------
    # struct -- Crystal structure object
    # max_wv (float) -- maximum wavevector
    # pw (float) -- peak width, units of A^-1
    # bg (float) -- average of normal background
    
    # returns lists of calcualted Q and intensity values
    #--------------------------------------------------------------------------
    
    q, intensity = struct.Scatter.generate_powder(max_wv, peak_width=pw, 
                                                  background=bg, 
                                                  powder_average=True)
    return q, intensity
    

#-----------------------------------------------------------------------------
''' normalize intensity values '''

def norm(ints):
    #--------------------------------------------------------------------------
    # ints -- list of calculated intensity values
    
    # returns list of normalilzed intensity values
    #--------------------------------------------------------------------------
    
    norm_ints = (ints / np.max(ints))
    
    return norm_ints



