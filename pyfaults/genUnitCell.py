##################################################################################
# Author: Sinclair R. Combs
##################################################################################

def genUnitCell(CSVPath, CSVName, lattParams, parentLyrNames,
                childLyrs=None):
    '''
    Parameters
    ----------
    CSVPath
        str : CSV file directory
    CSVName
        str : CSV file name
    lattParams
        nparray : lattice parameter values [a, b, c, alpha, beta, gamma]
    parentLyrNames
        list (str) : unique identifiers for each layer in unit cell
    childLyrs
        list (str, nparray) [optional] : information to generate child layers, 
        each list entry consists of the child layer name (str), name of the
        parent layer used to generate the child (str), and the translation
        vector applied to the parent (nparray)

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
    
    lyrs = pf.layer.getLayers(csv, latt, parentLyrNames)
    
    newLyrs = []
    if childLyrs is not None:
        for i in range(len(childLyrs)):
            childName = childLyrs[i][0]
            parentName = childLyrs[i][1]
            transVec = childLyrs[i][2]
            
            for j in range(len(lyrs)):
                if lyrs[j].layerName == parentName:
                    newLyr = lyrs[j].genChildLayer(childName, transVec)
                    newLyrs.append(newLyr)
                    
    for i in range(len(newLyrs)):
        lyrs.append(newLyrs[i])
        
    unitcell = pf.unitcell.Unitcell(CSVName, lyrs, latt)
    unitcell.toCif(CSVPath)
    
    return unitcell
        