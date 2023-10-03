##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import Dans_Diffraction as df
import copy as cp
import os

''' XRD SIMULATION METHODS '''
#---------------------------------------------------------------------------------
# simulates XRD pattern ----------------------------------------------------------
#---------------------------------------------------------------------------------
def fullSim(path, cif, wl, tt_max, pw=None, bg=None):
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

    Returns
    -------
    q : nparray
    ints : nparray
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
    return q, ints

#---------------------------------------------------------------------------------
# saves simulation as text file --------------------------------------------------
#---------------------------------------------------------------------------------
def saveSim(path, fn, q, ints):
    '''
    Parameters
    ----------
    path
        str : file save directory
    fn
        str : file save name
    q
        nparray : calculated Q values
    ints
        nparray : calculated intensity values
    '''
    with open(path + fn + '.txt', 'w') as f:
        for (q, ints) in zip(q, ints):
            f.write('{0} {1}\n'.format(q, ints))
    f.close()    
    
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

#---------------------------------------------------------------------------------
# calculates (hkl) reflections ---------------------------------------------------
#---------------------------------------------------------------------------------
def hklSim(path, cif, wl, tt_max):
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

    Returns
    -------
    reflections
        nparray : (hkl) reflections, 2theta, and intensity values
    '''
    struct = df.Crystal(path + cif + '.cif')
    
    # setup diffraction parameters
    e_kev = df.fc.wave2energy(wl)
    struct.Scatter.setup_scatter(scattering_type='xray', 
                                 energy_kev=e_kev, 
                                 min_twotheta=0, 
                                 max_twotheta=tt_max)
    
    # generate (hkl) reflections
    reflections = struct.Scatter.print_all_reflections(min_intensity=0.1)
    return reflections