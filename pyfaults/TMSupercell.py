##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import copy as cp
import numpy as np
import random as r

#---------------------------------------------------------------------------------
# TMSupercell object class
#---------------------------------------------------------------------------------
class TMSupercell(object):    
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
                 doc='fltLayer : faulted layer name')
    
    stackProb =\
        property(lambda self: self._stackProb, 
                 lambda self, val: self.setParam(stackProb=val),
                 doc='float : stacking probability')
        
    transMatrix =\
        property(lambda self: self._transitions,
                 lambda self, val: self.setParam(transitions=val))
    
    # creates instance of Supercell object ---------------------------------------
    def __init__(self, unitcell, nStacks, fltLayer, stackProb, transMatrix):
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
        self._stackProb = None
        self.setParam(nStacks)
        self.setLayers(unitcell, fltLayer, stackProb, transMatrix)
        return
    
    # set cell parameters --------------------------------------------------------
    def setParam(self, nStacks=None):
        if nStacks is not None:
            self._nStacks = nStacks
        return
    
    # set supercell layers --------------------------------------------------------
    def setLayers(self, nStacks, fltLayer, stackProb, transMatrix):
        
        seq = []
        seq.append(transMatrix['Start Layer'][0])
        
        numUCLyrs = 2

        for n in range(numUCLyrs * nStacks):
            p = (r.randint(0,100) * 0.01)

            possLyrIndex = []
            for i in range(len(transMatrix.index)):
                lyrName = transMatrix['Start Layer'][i]
                if lyrName == seq[-1]:
                    possLyrIndex.append(i)
        
            for l in possLyrIndex:
                isNext = False
            
                lyrProb = transMatrix['P'][l]

                if isinstance(lyrProb, str):
                    if lyrProb == 'P' or lyrProb == 'p':
                        if p <= stackProb:
                            isNext = True
                        elif '-' in lyrProb:
                            sub = lyrProb.split('-')
                            if p > stackProb and p < float(sub[0]):
                                isNext = True
                        elif float(lyrProb) == 1.0:
                            isNext = True
                        
                if isNext == True:
                    seq.append(transMatrix['Next Layer'][l])
                    break

        newLyrs = []
        newLyrs.append(lyrs[0])
    
        for i in range(1, len(seq)):
            lyrName = seq[i]
            if lyrName == 'F':
                lyrName = fltLayer

            if i%2 == 0:
                n = i/2
            elif i%2 != 0:
                n = (i+1)/2

            adj = []
            for j in range(len(transMatrix.index)):
                if transMatrix['Start Layer'][j] == lyrName:
                    if transMatrix['Next Layer'][j] == seq[i+1]:
                        x = transMatrix['x'][j]
                        y = transMatrix['y'][j]
                        z = transMatrix['z'][j]
                        adj = [x, y, z]
        
            newLyrAtoms = []
            for atom in lyrDict[lyrName]:
                for a in atomDict:
                    if a==atom:
                        newAtom = pf.layerAtom.LayerAtom(lyrName, atom, atomDict[a][0],
                                                         [float(atomDict[a][1]) + adj[0], 
                                                          float(atomDict[a][2]) + adj[1], 
                                                          float(atomDict[a][3]) + adj[2]], 
                                                         float(atomDict[a][4]), 
                                                         float(atomDict[a][5]), latt)
                        newLyrAtoms.append(newAtom)
            newLyr = pf.layer.Layer(newLyrAtoms, latt, lyrName)
            newLyrs.append(newLyr)
            
        self._layers = newLyrs
        return
    
    

