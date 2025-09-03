##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import Dans_Diffraction as df

# simulates XRD pattern, returns normalized intensity values ----------
def fullSim(path, cif, wl, tt_max, savePath, pw=None, bg=None):
    '''
    Parameters
    ----------
    path (str) : CIF file location
    cif (str) : CIF file name
    wl (float) : instrument wavelength
    tt_max (float) : maximum two theta
    savePath (str) : location to save simulation data to
    pw (float, optional) : peak broadening term
    bg (float, optional) : average of normal background

    Returns
    -------
    q (array_like) : calculated Q values
    ints (array_like) : calculated intensity values, normalized
    '''
    
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
    
    # save simulation data
    with open(savePath + cif + '_sim.txt', 'w') as f:
        for (q, ints) in zip(q, ints):
            f.write('{0} {1}\n'.format(q, ints))
    f.close() 
    
    return q, ints 
