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

    Returns
    -------
    simData
        DataFrame : tabulated simulation data, includes the following as columns
            'Model' -- unique identifier
            'Sim Q' -- simulated Q values
            'Sim Norm Intensity' -- simulated intensity values, normalized
    '''
    import pyfaults as pf
    import os
    import pandas as pd
    
    # create 'sims' folder in file directory
    if os.path.exists(savePath + 'sims/') == False:
        os.mkdir(savePath + 'sims/')
        
    df = pd.read_csv(CIFPath + CSVName + '.csv')
    
    simQList = []
    simDiffList =[]
    
    # simulate XRD pattern for each supercell in directory
    for i in df.index:
        name = df['Model'][i]
        q, ints = pf.simXRD.fullSim(CIFPath, name, wl, maxTT, pw=pw)
        simQList.append(q)
        simDiffList.append(pf.norm(ints))
        # save simulation data
        with open(savePath + 'sims/' + name + '_sim.txt', 'w') as f:
            for i in range(len(q)):
                f.write('{0} {1}\n'.format(q[i], ints[i]))
        f.close()
    
    # store all simulations in data frame
    simData = pd.dataFrame()
    simData['Model'] = df['Model']
    simData['Sim Q'] = simQList
    simData['Sim Norm Intensity'] = simDiffList
    
    simData.to_csv(savePath + 'sims/sim_info.csv')
    
    return simData