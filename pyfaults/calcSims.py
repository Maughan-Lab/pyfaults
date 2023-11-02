##################################################################################
# Author: Sinclair R. Combs
##################################################################################

def calcSims(CIFPath, CSVName, wl, maxTT, pw, savePath):
    '''
    Parameters
    ----------
    CIFPath
        str : supercell file directory
    CSVName
        str : file name of metadata CSV generated from genSupercells method
    wl
        float : instrument wavelength (A)
    maxTT
        float : maximum 2theta value (degrees)
    pw
        float : simulated peak broadening term
    savePath
        str : file directory to save 'sims' folder
    '''
    import pyfaults as pf
    import os
    import pandas as pd
    
    # create 'sims' folder in file directory
    if os.path.exists(savePath + 'sims/') == False:
        os.mkdir(savePath + 'sims/')
        
    df = pd.read_csv(CIFPath + CSVName + '.csv')
    
    # simulate XRD pattern for each supercell in directory
    for i in df.index:
        name = df['Model'][i]
        q, ints = pf.simXRD.fullSim(CIFPath, name, wl, maxTT, pw=pw, 
                                    savePath=savePath + 'sims/')
    
    return
