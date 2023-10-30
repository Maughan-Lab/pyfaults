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
def fullSim(path, cif, wl, tt_max, pw=None, bg=None, norm=False):
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
    pw
        float (optional) : peak width (A^-1)
    bg
        float (optional) : average of normal background
    norm
        bool (optional) : set to True if intensities are to be normalized
        The d

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
    with open(path + cif + '_sim.txt', 'w') as f:
        for (q, norm_ints) in zip(q, norm_ints):
            f.write('{0} {1}\n'.format(q, norm_ints))
    f.close() 
    
    return q, norm_ints  
    
#---------------------------------------------------------------------------------
# calculates simulated diffraction pattern for DataFrame of supercells -----------
#---------------------------------------------------------------------------------
def calcSims(df, path, wl, maxTT, pw):
    '''
    Parameters
    ----------
    df
        DataFrame : tabulated data returned from genSupercells()
    path
        str : file save directory
    wl
        float : instrument wavelength (A)
    maxTT
        float : maximum 2theta value (degrees)
    pw
        float : simulated peak broadening term

    Returns
    -------
    simDF
        DataFrame : tabulated simulation data, includes the following
            'Model' -- unique identifier
            'Stacking Vector' -- stacking vector in string format
            'Stacking Probability' -- stacking fault probability in str format
            'S_x' -- x-component of stacking vector
            'S_y' -- y-component of stacking vector
            'S_z' -- z-component of stacking vector
            'P' -- fault probability (0 to 1)
            'Simulated Q' -- calculated Q values (A^-1)
            'Simulated Intensity' -- calculated intensity values
    '''
    from pyfaults import norm
    # create 'sims' folder in file directory
    if os.path.exists(path + 'sims/') == False:
        os.mkdir(path + 'sims/')
    
    # create copy of simDF dataframe
    simDF = cp.deepcopy(df)
    
    simQList = []
    simDiffList = []
    # simulate XRD pattern for each entry in simDF
    for i in simDF.index:
        name = simDF['Model'][i]
        Q, ints = fullSim(path, name, wl, maxTT, pw=pw)
        simQList.append(Q)
        simDiffList.append(norm(ints))
        # save simulation data
        with open(path + 'sims/' + name + '_sim.txt', 'w') as f:
            for (Q, ints) in zip(Q, ints):
                f.write('{0} {1}\n'.format(Q, ints))
        f.close()
     
    # append each entry with XRD simulation data
    simDF['Simulated Q'] = simQList
    simDF['Simulated Intensity'] = simDiffList
    return simDF