##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import Dans_Diffraction as df
import copy as cp
import os

''' XRD SIMULATION METHODS '''
#---------------------------------------------------------------------------------
# simulates XRD pattern, returns normalized intensity values ---------------------
#---------------------------------------------------------------------------------
def fullSim(path, cif, wl, tt_max, savePath, pw=None, bg=None):
    '''
    Parameters
    ----------
    path
        str : CIF file directory
    cif
        str : CIF file name
    wl
        float : instrument wavelength (A)
    tt_max
        float : maximum 2theta (degrees)
    savePath
        str : file directory to save simulation data
    pw
        float (optional) : peak width (A^-1)
    bg
        float (optional) : average of normal background

    Returns
    -------
    q : nparray
    ints : nparray
    '''
    from pyfaults import norm
    
    # load cif
    struct = df.Crystal(path + cif + '.cif')
    
    # setup diffraction parameters
    energy_kev = df.fc.wave2energy(wl)
    struct.Scatter.setup_scatter('xray')
    wavevector_max = df.fc.calqmag(tt_max, energy_kev)
    if pw is None:
        pw = 0.0
    if bg is None:
        bg = 0
    
    # simultate XRD pattern
    q, ints = struct.Scatter.generate_powder(wavevector_max, 
                                             peak_width=pw, 
                                             background=bg, 
                                             powder_average=True)
    
    norm_ints = norm(ints)
    
    # save simulation data
    with open(savePath + cif + '_sim.txt', 'w') as f:
        for (q, norm_ints) in zip(q, norm_ints):
            f.write('{0} {1}\n'.format(q, norm_ints))
    f.close() 
    
    return q, norm_ints 