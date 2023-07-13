"""
XRD Simulations
"""

import os
import Dans_Diffraction as df
import numpy as np

#-----------------------------------------------------------------------------
def full_XRD_sim(cif_dir, cif, wavelength, tt_max, pw, bg):
    # cif_dir (str) -- path to folder where cif is stored
    # cif (str) -- name of cif file
    # wavelength (float) -- wavelength of instrument used to collect data
    # tt_max (float) -- maximum 2theta value for simulation
    # pw (float) -- peak width, in A^-1
    # bg (float) -- average of normal background
    
    #load cif
    path = os.path.join(cif_dir, cif + ".cif")
    struct = df.Crystal(path)
    
    # setup
    energy_kev = df.fc.wave2energy(wavelength)
    struct.Scatter.setup_scatter("xray")
    wavevector_max = df.fc.calqmag(tt_max, energy_kev)
    
    # simultate XRD
    q, ints = struct.Scatter.generate_powder(wavevector_max, 
                                             peak_width=pw, background=bg,
                                             powder_average=True)
    return q, ints
    # returns lists of calcualted Q and intensity values


#-----------------------------------------------------------------------------
# save .txt file of simulated diffraction information
def save_sim(directory, name, q, ints):
    # directory (str) -- path to folder where file is to be saved
    # name (str) -- name of generated .txt file
    # q (list) -- calculated Q value list
    # ints (list) -- calculated intensity value list
    
    path = os.path.join(directory, name + ".txt")
    
    new_file = open(path, "w")
    
    with new_file as f:
        for (q, ints) in zip(q, ints):
            f.write("{0} {1}\n".format(q, ints))
    
    new_file.close()
    
    
#-----------------------------------------------------------------------------
# import XRD simulation .txt file
def import_sim(path, name):
    # path (str) -- path to folder where text file is stored
    # name (str) -- name of text file
    
    filepath = path + name + ".txt"
    q, ints = np.loadtxt(filepath, unpack=True, dtype=float)
    
    return q, ints
    # returns lists of Q and intensity values from text file
    

#-----------------------------------------------------------------------------
# load cif as structure
def load_cif(directory, filename):
    # directory (str) -- path to folder where cif is stored
    # filename (str) -- name of cif file
    
    path = os.path.join(directory, filename + ".cif")
    return df.Crystal(path)
    # returns structure as Crystal object


#-----------------------------------------------------------------------------
# set up diffraction simulation conditions
def sim_setup(struct, wl, max_tt):
    # struct (Crystal) -- crystal structure to be simulated
    # wl (float) -- wavelength of instrument used to collect data
    # max_tt (float) -- maximum 2theta value for simulation
    
    energy_kev = df.fc.wave2energy(wl)
    struct.Scatter.setup_scatter("xray")
    max_wavevector = df.fc.calqmag(max_tt, energy_kev)
    return max_wavevector
    # returns maximum wavevector
    

#-----------------------------------------------------------------------------
# simulate powder diffraction
def sim(struct, max_wv, pw, bg):
    # struct (Crystal) -- crystal structure to be simulated
    # max_wv (float) -- maximum wavevector, calculated in sim_setup
    # pw (float) -- peak width, in A^-1
    # bg (float) -- average of normal background
    
    q, intensity = struct.Scatter.generate_powder(max_wv, peak_width=pw, 
                                                  background=bg, 
                                                  powder_average=True)
    return q, intensity
    # returns lists of calcualted Q and intensity values
    

#-----------------------------------------------------------------------------
# normalize intensity values
def norm(ints):
    # ints (list) -- calculated intensity value list
    
    norm_ints = (ints / np.max(ints))
    return norm_ints
    # returns list of normalized intensity values



