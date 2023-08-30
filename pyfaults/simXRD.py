#########################################################################################
# pyfaults.simXRD
# Author: Sinclair R. Combs
#########################################################################################

import Dans_Diffraction as df
import numpy as np

#########################################################################################
''' XRD simulation functions '''
#########################################################################################

def fullSim(path, cif, wl, tt_max, pw=None, bg=None):
    '''
    Simulates XRD pattern

    Parameters
    ----------
    path : str
        File directory where cif is store
    cif : str
        File name
    wl : float
        Instrument wavelength (A)
    tt_max : float
        Maximum 2theta (deg)
    pw : float, optional
        Peak width (A^-1)
    bg : float, optional
        Average of normal background

    Returns
    -------
    q : nparray
        Calculated Q values
    ints : nparray
        Calculated intensity values
    '''
    
    # load cif
    struct = df.Crystal(path + cif + ".cif")
    
    # setup diffraction parameters
    energy_kev = df.fc.wave2energy(wl)
    struct.Scatter.setup_scatter("xray")
    wavevector_max = df.fc.calqmag(tt_max, energy_kev)
    
    if pw is None:
        pw = 0.0
    if bg is None:
        bg = 0
    
    # simultate XRD pattern
    q, ints = struct.Scatter.generate_powder(wavevector_max, peak_width=pw, 
                                             background=bg, powder_average=True)
    return q, ints


#----------------------------------------------------------------------------------------
def saveSim(path, fn, q, ints):
    '''
    Saves simulated XRD as text file

    Parameters
    ----------
    path : str
        File directory
    fn : str
        File name
    q : nparray
        Calculated Q values
    ints : nparray
        Calculated intensity values
    '''
    with open(path + fn + ".txt", "w") as f:
        for (q, ints) in zip(q, ints):
            f.write("{0} {1}\n".format(q, ints))
    f.close()    

#----------------------------------------------------------------------------------------
def loadCif(path, fn):
    '''
    Imports cif file as diffpy Structure object

    Parameters
    ----------
    path : str
        File directory
    fn : str
        File name

    Returns
    -------
    struct : Structure
        Structure object

    '''
    struct = df.Crystal(path + fn + ".cif")
    return struct

#----------------------------------------------------------------------------------------
def norm(ints):
    '''
    Normalize intensity values

    Parameters
    ----------
    ints : nparray
        Intensity values

    Returns
    -------
    norm_ints : nparray
        Normalized intensity values
    '''
    norm_ints = (ints / np.max(ints))
    return norm_ints


#----------------------------------------------------------------------------------------
def hklSim(path, cif, wl, tt_max):
    '''
    Calculates (hkl) reflections

    Parameters
    ----------
    path : str
        File directory where cif is store
    cif : str
        File name
    wl : float
        Instrument wavelength (A)
    tt_max : float
        Maximum 2theta (deg)

    Returns
    -------
    reflections : nparray
        Calculated (hkl) reflections, 2theta, and intensity values
    '''
    # load cif
    struct = df.Crystal(path + cif + ".cif")
    
    # setup diffraction parameters
    e_kev = df.fc.wave2energy(wl)
    struct.Scatter.setup_scatter(scattering_type="xray", energy_kev=e_kev, 
                                 min_twotheta=0, max_twotheta=tt_max)
    
    # generate (hkl) reflections
    reflections = struct.Scatter.print_all_reflections(min_intensity=0.1)
    
    return reflections


#----------------------------------------------------------------------------------------
def savehkl(path, fn, data):
    '''
    Saves (hkl) data as text file

    Parameters
    ----------
    path : str
        File directory
    fn : str
        File name
    data : nparray
        (hkl) data
    '''
    with open(path + fn + "_hkl.txt", "w") as f:
        f.write(data)
    f.close()
    
#----------------------------------------------------------------------------------------
def diffCurve(q1, q2, ints1, ints2): 
    '''
    Calculates difference curve between two datasets

    Parameters
    ----------
    q1 : nparray
        Q data from dataset one
    q2 : nparray
        Q data from dataset two
    ints1 : nparray
        Intensity data from dataset one
    ints2 : nparray
        Intensity data from dataset two

    Returns
    -------
    diff_q : noarray
        Q data of difference curve
    diff_ints : nparray
        Intensity data of difference curve

    '''
    diff_q_list = []
    diff_ints_list = []
    for i in range(0, len(q1)):
        q1_val = float('%.3f'%(q1[i]))
        for j in range(0, len(q2)):
            q2_val = float('%.3f'%(q2[j]))
            if q1_val == q2_val:
                diff_q_list.append(q1_val)
                diff_ints_list.append(ints1[i]-ints2[j])
    diff_q = np.array(diff_q_list)
    diff_ints = np.array(diff_ints_list)
    
    return diff_q, diff_ints


#----------------------------------------------------------------------------------------
def saveDiffCurve(q, ints, path, fn):
    '''
    Saves difference curve to text file

    Parameters
    ----------
    q : nparray
        Q data of difference curve
    ints : nparray
        Intensity data of difference curve
    path : str
        File directory
    fn : str
        File name
    '''
    with open(path + fn + ".txt", "w") as f:
        for (q, ints) in zip(q, ints):
            f.write("{0} {1}\n".format(q, ints))
    f.close()