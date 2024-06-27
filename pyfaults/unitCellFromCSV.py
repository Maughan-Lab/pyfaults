##################################################################################
# Author: Sinclair R. Combs
##################################################################################

#---------------------------------------------------------------------------------
# generates unit cell from CSV of atomic parameters ------------------------------
#---------------------------------------------------------------------------------
def unitCellFromCSV(CSVPath, CSVName, lattParams, lyrNames):
    '''
    Parameters
    ----------
    CSVPath
        str : CSV file directory
    CSVName
        str : CSV file name
    lattParams
        nparray : lattice parameter values [a, b, c, alpha, beta, gamma]
    lyrNames
        list (str) : unique identifiers for each layer in unit cell

    Returns
    -------
    unitcell
        Unitcell : layered unit cell generated from CSV parameters
    '''
    import pyfaults as pf
    
    csv = pf.importCSV(CSVPath, CSVName)
    
    latt = pf.lattice.Lattice(a=lattParams[0],
                              b=lattParams[1],
                              c=lattParams[2],
                              alpha=lattParams[3],
                              beta=lattParams[4],
                              gamma=lattParams[5])
    
    lyrs = pf.layer.getLayers(csv, latt, lyrNames)
    
    unitcell = pf.unitcell.Unitcell(CSVName, lyrs, latt)
    unitcell.toCif(CSVPath)
    
    return unitcell
        