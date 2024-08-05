##################################################################################
# Author: Sinclair R. Combs
##################################################################################

#---------------------------------------------------------------------------------
# generates unit cell from CSV of atomic parameters ------------------------------
#---------------------------------------------------------------------------------
def importCSV(path, fn, lattParams, lyrNames):
    
    import pyfaults as pf

    csv = pd.read_csv(path + fn + '.csv')
    
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
