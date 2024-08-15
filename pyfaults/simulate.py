import os, glob

def simulateRev(path):
    '''
    Parameters
    ----------
    path (str) : input file path
    '''
    
    import pyfaults as pf
    
    unitcell, ucDF, gsDF, scDF, simDF = pf.pfInput.pfInput(path)

    wl = simDF.loc[0, 'wl']
    maxTT = simDF.loc[0, 'maxTT']
    pw = simDF.loc[0, 'pw']
    
    # creates folder to store simulation data
    if os.path.exists('./simulations') == False:
        os.mkdir('./simulations')

    fileList = glob.glob('./supercells/*')
    for f in range(len(fileList)):
        fileList[f] = fileList[f].replace('.cif', '')
        fileList[f] = fileList[f].replace('./supercells\\', '')
    
    # PXRD simulation for each supercell
    for f in fileList:
        q, ints = pf.simXRD.fullSim('./supercells/', f, wl.iloc[0], maxTT.iloc[0], pw=pw.iloc[0], savePath='./simulations/')
        