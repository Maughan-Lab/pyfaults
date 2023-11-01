##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import copy as cp
import numpy as np
import random as r
import pandas as pd
import os

''' SUPERCELL CLASS '''
#---------------------------------------------------------------------------------
class Supercell(object):    
    # properties -----------------------------------------------------------------
    unitcell =\
        property(lambda self: self._unitcell,
                 doc='Unitcell : base unit cell')
    
    lattice =\
        property(lambda self: self._lattice,
                 doc='Lattice : unit cell lattice parameters')
    
    nStacks =\
        property(lambda self: self._nStacks, 
                 lambda self, val: self.setParam(nStacks=val),
                 doc='int : number of stacks in supercell')
    
    layers =\
        property(lambda self: self._layers, 
                 lambda self, val: self.setLayers(layers=val),
                 doc='list (Layer) : layers contained in supercell')
    
    fltLayer =\
        property(lambda self: self._fltLayer,
                 lambda self, val: self.setParam(fltLayer=val),
                 doc='fltLayer [optional] : faulted layer name')
    
    stackVec =\
        property(lambda self: self._stackVec, 
                 lambda self, val: self.setParam(stackVec=val),
                 doc='nparray [optional] : stacking vector')
    
    stackProb =\
        property(lambda self: self._stackProb, 
                 lambda self, val: self.setParam(stackProb=val),
                 doc='float [optional] : stacking probability')
    
    # creates instance of Supercell object ---------------------------------------
    def __init__(self, unitcell, nStacks, 
                 fltLayer=None, stackVec=None, stackProb=None):
        from pyfaults.lattice import Lattice
        # initialize parameters
        self._unitcell = unitcell
        newLatt = Lattice(unitcell.lattice.a,
                          unitcell.lattice.b,
                          (unitcell.lattice.c * nStacks),
                          unitcell.lattice.alpha,
                          unitcell.lattice.beta,
                          unitcell.lattice.gamma)
        self._lattice = newLatt
        self._nStacks = None
        self._layers = None
        self._fltLayer = None
        self._stackVec = None
        self._stackProb = None
        self.setParam(nStacks)
        self.setLayers(unitcell, fltLayer, stackVec, stackProb)
        return
    
    # set cell parameters --------------------------------------------------------
    def setParam(self, nStacks=None):
        if nStacks is not None:
            self._nStacks = nStacks
        return
    
    # set supercell layers --------------------------------------------------------
    def setLayers(self, unitcell, fltLayer=None, stackVec=None, stackProb=None):
        newLayers = []
        for lyr in self.unitcell.layers:
            for n in range(self.nStacks):
                tag = '_n' + str(n+1)
                p = r.randint(0,100)
                
                if fltLayer is not None:
                    self._fltLayer = fltLayer
                    self._stackVec = stackVec
                    self._stackProb = stackProb
                    if lyr.layerName == fltLayer:
                        if p <= (stackProb * 100):
                            isFlt = True
                        elif p > (stackProb * 100):
                            isFlt = False
                    elif lyr.layerName != fltLayer:
                        isFlt = False
                elif fltLayer is None: 
                    isFlt = False
                
                newLyr = cp.deepcopy(lyr)
                
                if isFlt == True:
                    newLayerName = lyr.layerName + tag + '_fault'
                    newLyr.setParam(layerName=newLayerName, 
                                    lattice=self.lattice)
                    for atom in newLyr.atoms:
                        alabel = atom.atomLabel.split('_')
                        newXYZ = [atom.x, 
                                  atom.y, 
                                  ((atom.z + n) / self.nStacks)]
                        fltXYZ = np.add(newXYZ, stackVec)
                        atom.setParam(layerName=newLayerName, 
                                      atomLabel=alabel[0], 
                                      xyz=fltXYZ, 
                                      lattice=self.lattice)
                if isFlt == False:
                    newLayerName = lyr.layerName + tag
                    newLyr.setParam(layerName=newLayerName, 
                                    lattice=self.lattice)
                    for atom in newLyr.atoms:
                        alabel = atom.atomLabel.split('_')
                        newXYZ = [atom.x, 
                                  atom.y, 
                                  ((atom.z + n) / self.nStacks)]
                        atom.setParam(layerName=newLayerName, 
                                      atomLabel=alabel[0], 
                                      xyz=newXYZ, 
                                      lattice=self.lattice)
                
                newLayers.append(newLyr)
        self._layers = newLayers
        return
    
    # prints names of faulted layers ---------------------------------------------
    def show_faults(self):
        for lyr in self.layers:
            if "fault" in lyr.layerName:
                print(lyr.layerName)
        return
                
    # generates CIF of layer -----------------------------------------------------
    def toCif(self, path):
        from pyfaults import toCif
        name = self.unitcell.name
        toCif(self, path, 'Supercell_' + name + '_N' + str(self.nStacks) + '_CIF')
                
''' SUPERCELL GENERATOR METHOD '''                
#---------------------------------------------------------------------------------
# generates supercells within search parameter grid ------------------------------
#---------------------------------------------------------------------------------
def genSupercells(unitcell, nStacks, fltLayer, probList, sVecList, path):
    '''
    Parameters
    ----------
    unitcell
        Unitcell : base unit cell stack
    nStacks
        int : number of stacks per supercell
    fltLayer
        str : name of faulted layer
    probList
        list (float) : stacking fault probabilities to search
    sVecList
        list (nparray) : stacking vectors to search
    path
        str : CIF file save directory

    Returns
    -------
    df
        DataFrame : tabulated supercell data, includes the following as columns
            'Model' -- unique identifier
            'Stacking Vector' -- stacking vector in string format
            'Stacking Probability' -- stacking fault probability in str format
            'S_x' -- x-component of stacking vector
            'S_y' -- y-component of stacking vector
            'S_z' -- z-component of stacking vector
            'P' -- fault probability (0 to 1)
    '''
    from pyfaults import toCif, formatStr
    # create 'CIFs' folder in file directory
    if os.path.exists(path + 'CIFs/') == False:
        os.mkdir(path + 'CIFs/')
    
    df = pd.DataFrame()
    cellList = []
    probStrList, sVecStrList = formatStr(probList, sVecList)
    
    # generate unfaulted supercell
    UF = Supercell(unitcell, nStacks)
    cellList.append([UF, 'Unfaulted'])
    modelCol = ['Unfaulted']
    sVecTxtCol = [np.nan]
    probTxtCol = [np.nan]
    sxCol = [np.nan]
    syCol = [np.nan]
    szCol = [np.nan]
    probCol = [np.nan]
    
    # generate each faulted supercell
    for p in range(len(probList)):
        for s in range(len(sVecList)):
            FLT = Supercell(unitcell, nStacks, fltLayer=fltLayer, 
                            stackVec=sVecList[s], stackProb=probList[p])
            cellTag = 'S' + str(s+1) + '_P' + str(int(probList[p]*100))
            cellList.append([FLT, cellTag])
            modelCol.append(cellTag)
            sVecTxtCol.append(sVecStrList[s])
            probTxtCol.append(probStrList[p])
            sxCol.append(sVecList[s][0])
            syCol.append(sVecList[s][1])
            szCol.append(sVecList[s][2])
            probCol.append(probList[p])
    
    # export CIF for each supercell
    for c in range(len(cellList)):
        toCif((cellList[c][0]), path + 'CIFs/', cellList[c][1])
    
    # add entry for each supercell
    df['Model'] = modelCol
    df['Stacking Vector'] = sVecTxtCol
    df['Stacking Probability'] = probTxtCol
    df['S_x'] = sxCol
    df['S_y'] = syCol
    df['S_z'] = szCol
    df['P'] = probCol
    return df               