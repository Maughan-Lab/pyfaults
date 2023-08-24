'''
UF -- collection of Unitcell objects constructing an unfaulted supercell
FLT -- collection of Unitcell objects constructing a faulted supercell
'''

import copy as cp
import numpy as np
import random as r

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' UF CLASS '''
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class UF(list):
    def __init__(self, unitcell, n_stacks):
        '''
        UF class initialization

        Parameters
        ----------
        unitcell : Unitcell
            Single unit cell stack
        n_stacks : int
            Number of stacks in supercell

        Returns
        -------
        None

        '''

        # assign n_stacks variable
        self.n_stacks = n_stacks
        
        # create a copy of the unit cell lattice
        latt = cp.deepcopy(unitcell.lattice)
        
        # scale the c-axis vector by the number of stacks
        latt.setLatPar(c = n_stacks * latt.c)
        
        # create and assign lattice variable
        self.lattice = latt
        
        # create a copy of the unit cell for each N in n_stacks
        for n in range(0, n_stacks):
            uc = cp.deepcopy(unitcell)
            
            # create a copy of each layer in the unit cell
            for layer in uc:
                wl = cp.deepcopy(layer)
                
                # scale the z-position of each atom and add label based on N
                for a in wl:
                    a.z = (a.z + n) / n_stacks
                    a.label = a.label + "_n" + str(n)
                
                # set lattice of layer to N-scaled lattice
                wl.lattice = latt
                
                # add label to layer name based on N
                wl.layer_name = wl.layer_name + "_n" + str(n)
                
                # append layers to the UF_supercell object (list)
                self.append(wl)

        return
    
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    ''' UF CLASS METHODS '''
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    
    def show_atoms(self):
        '''
        Prints atom labels and [x,y,z] of all atoms in layer

        Returns
        -------
        None

        '''
        for l in self:
            for a in l:
                print(a.label, a.xyz)
    
    #--------------------------------------------------------------------------
    def show_layers(self):
        '''
        Prints labels all layers in unit cell

        Returns
        -------
        None

        '''
        layer = self.layers
        for i in layer:
            i.layer_name


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' FLT CLASS '''
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class FLT_supercell(list):
    def __init__(self, unitcell, n_stacks, flt_layer, shift, p):
        '''
        FLT class initialization

        Parameters
        ----------
        unitcell : Unitcell
            Single unit cell stack
        n_stacks : int
            Number of stacks in supercell
        flt_layer : str
            Name of faulted layer
        shift : array
            Stacking vector [x, y, 0]
        p : float
            Stacking probability

        Returns
        -------
        None

        '''

        # assign n_stacks, shift, and p variables
        self.n_stacks = n_stacks
        self.shift = shift
        self.p = p
        
        # create a copy of the unit cell lattice
        latt = cp.deepcopy(unitcell.lattice)
        
        # scale the c-axis vector by the number of stacks
        latt.setLatPar(c = n_stacks * latt.c)
        
        # create and assign lattice variable
        self.lattice = latt
        
        # convert decimal probability to a percentage
        percent = p * 100
        
        # create a copy of the unit cell for each N in n_stacks
        for n in range(0, n_stacks):
            n_label = "_n" + str(n)
            flt_cell = cp.deepcopy(unitcell)
            
            for layer in flt_cell:
                # generates a random integer to assess probability of fault
                gen_prob = r.randint(0, 100)
                
                # if layer matches the name of the given faulted layer and the
                # probability condition is met
                if layer.layer_name == flt_layer and gen_prob <= percent:
                    
                    # create a copy of faulted layer
                    flt = cp.deepcopy(layer)
                    
                    # add stacking vector, scale z-position, and add labels to
                    # each atom in the faulted layer
                    for a in flt:
                        a.xyz = np.add(a.xyz, shift)
                        a.z = (a.z + n) / n_stacks
                        a.label = a.label + "_flt" + n_label
                        
                    # set lattice of layer to N-scaled lattice
                    flt.lattice = latt
                    
                    # add label to layer name based on N
                    flt.layer_name = flt.layer_name + "_flt" + n_label
                    
                    # append faulted layer to the FLT_supercell object (list)
                    self.append(flt)
                    
                # if layer matches the name of the given faulted layer and the
                # probability condition is NOT met
                elif layer.layer_name == flt_layer and gen_prob > percent:
                    
                    # create a copy of layer
                    flt = cp.deepcopy(layer)
                    
                    # scale z-position and add labels to each atom in layer
                    for a in flt:
                        a.z = (a.z + n) / n_stacks
                        a.label = a.label + n_label
                        
                    # set lattice of layer to N-scaled lattice
                    flt.lattice = latt
                    
                    # add label to layer name based on N
                    flt.layer_name = flt.layer_name + n_label
                    
                    # append layer to the FLT_supercell object (list)
                    self.append(flt)
                    
                # if layer does not match the name of the given faulted layer
                elif layer.layer_name != flt_layer:
                    
                    # create a copy of layer
                    l = cp.deepcopy(layer)
                    
                    # scale z-position and add labels to each atom in layer
                    for a in l:
                        a.z = (a.z + n) / n_stacks
                        a.label = a.label + n_label
                        
                    # set lattice of layer to N-scaled lattice
                    l.lattice = latt
                    
                    # add label to layer name based on N
                    l.layer_name = l.layer_name + n_label
                    
                    # append layer to the FLT_supercell object (list)
                    self.append(l)

        return
    
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    ''' FLT CLASS METHODS '''
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------

    def show_atoms(self):
        '''
        Prints atom labels and [x,y,z] of all atoms in layer

        Returns
        -------
        None

        '''
        for l in self:
            for a in l:
                print(a.label, a.xyz)
    
    #--------------------------------------------------------------------------
    def show_layers(self):
        '''
        Prints labels all layers in unit cell

        Returns
        -------
        None

        '''
        layer = self.layers
        for i in layer:
            i.layer_name


