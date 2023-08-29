#########################################################################################
# pyfaults.constructCell
# Author: Sinclair R. Combs
#########################################################################################

import copy as cp
import numpy as np
import pandas as pd
import random as r

#########################################################################################
''' Lattice object class -- Stores unit cell lattice parameters '''
#########################################################################################

class Lattice(object):
    '''
    Parameters
    ----------
    a : float
        Unit cell lattice vector a (A)
    b : float
        Unit cell lattice vector b (A)
    c : float
        Unit cell lattice vector c (A)
    alpha : float
        Unit cell lattice angle alpha (deg)
    beta : float
        Unit cell lattice angle beta (deg)
    gamma : float
        Unit cell lattice angle gamma (deg)
    '''
    
    # properties ------------------------------------------------------------------------
    a = property(lambda self: self._a,
                 lambda self, val: self.setParam(a=val),
                 doc="float : unit cell vector a")
    
    b = property(lambda self: self._b,
                 lambda self, val: self.setParam(b=val),
                 doc="float : unit cell vector b")
    
    c = property(lambda self: self._c,
                 lambda self, val: self.setParam(c=val),
                 doc="float : unit cell vector c")
    
    alpha = property(lambda self: self._alpha,
                 lambda self, val: self.setParam(alpha=val),
                 doc="float : unit cell angle alpha")
    
    beta = property(lambda self: self._beta,
                 lambda self, val: self.setParam(beta=val),
                 doc="float : unit cell angle beta")
    
    gamma = property(lambda self: self._gamma,
                 lambda self, val: self.setParam(gamma=val),
                 doc="float : unit cell angle gamma")
    
    # creates instance of Lattice object ------------------------------------------------
    def __init__(self, a, b, c, alpha, beta, gamma):
        
        # initialize parameters
        self._a = None
        self._b = None
        self._c = None
        self._alpha = None
        self._beta = None
        self._gamma = None
        
        self.setParam(a, b, c, alpha, beta, gamma)
        
        return
    
    # sets lattice parameters -----------------------------------------------------------
    def setParam(self, a=None, b=None, c=None, alpha=None, beta=None, gamma=None):
        
        if a is not None:
            self._a = a
        if b is not None:
            self._b = b
        if c is not None:
            self._c = c
        if alpha is not None:
            self._alpha = alpha
        if beta is not None:
            self._beta = beta
        if gamma is not None:
            self._gamma = gamma
            
        return
    
    # prints lattice parameters ---------------------------------------------------------
    def display(self):
        latt = "\n".join(("a : " + self.a, 
                       "b : " + self.b,
                       "c : " + self.c,
                       "alpha : " + self.alpha,
                       "beta : " + self.beta,
                       "gamma : " + self.gamma))
        print(latt)
        
    
#########################################################################################
''' LayerAtom object class -- Stores properties of atom in specific layer '''
#########################################################################################

class LayerAtom(object):
    '''
    Parameters
    ----------
    layerName : str
        Unique identifier for layer in which atom resides
    atomLabel : str
        Unique identifier for LayerAtom object
    element : str
        Atomic element and oxidation state if applicable
    xyz : nparray
        Atomic position in fractional coordinates
    occupancy : float
        Crystallographic site occupancy
    lattice : Lattice
        Unit cell lattice parameters
    '''
        
    # properties ------------------------------------------------------------------------
    layerName = property(lambda self: self._layerName,
                         lambda self, val: self.setParam(layerName=val))
        
    atomLabel = property(lambda self: self._atomLabel,
                         lambda self, val: self.setParam(atomLabel=val))
        
    element = property(lambda self: self._element,
                       lambda self, val: self.setParam(element=val))
        
    xyz = property(lambda self: self._xyz,
                   lambda self, val: self.setParam(xyz=val))
        
    x = property(lambda self: self._xyz[0],
                 lambda self, val: self.setParam(0, val),
                 doc="float : fractional coordinate x")
        
    y = property(lambda self: self._xyz[1],
                 lambda self, val: self.setParam(1, val),
                 doc="float : fractional coordinate y")
        
    z = property(lambda self: self._xyz[2],
                 lambda self, val: self.setParam(2, val),
                 doc="float : fractional coordinate z")
        
    occupancy = property(lambda self: self._occupancy,
                         lambda self, val: self.setParam(occupancy=val))
        
    lattice = property(lambda self: self._lattice,
                       lambda self, val: self.setParam(lattice=val))
    
    # creates instance of LayerAtom object ----------------------------------------------
    def __init__(self, layerName, atomLabel, element, xyz, occupancy, lattice):
        
        # initialize parameters
        self._layerName = None
        self._atomLabel = None
        self._element = None
        self._xyz = None
        self._x = None
        self._y = None
        self._z = None
        self._lattice = None
        self._occupancy = None
        
        self.setParam(layerName, atomLabel, element, xyz, lattice, occupancy)
        
        return
    
    # sets atom parameters --------------------------------------------------------------
    def setParam(self, layerName=None, atomLabel=None, element=None, xyz=None,
                 lattice=None, occupancy=None):
        
        if layerName is not None:
            self._layerName = layerName
        if atomLabel is not None:
            self._atomLabel = atomLabel + "_" + layerName
        if element is not None:
            self._element = element
        if xyz is not None:
            self._xyz = xyz
            self._x = xyz[0]
            self._y = xyz[1]
            self._z = xyz[2]
        if lattice is not None:
            self._lattice = lattice
        if occupancy is not None:
            self._occupancy = occupancy
        
        return
    
    # prints atom layer, name, element, position, and occupancy -------------------------
    def display(self):
        a = "\n".join(("Layer Name: " + self.layerName, 
                       "Atom Label: " + self.atomLabel,
                       "Element: " + self.element,
                       "Atomic Position: " + str(self.xyz),
                       "Occupancy: " + str(self.occupancy)))
        print(a)


#########################################################################################
''' Layer object class -- Stores layer properties'''
#########################################################################################

class Layer(object):
    '''
    Parameters
    ----------
    atoms : list (LayerAtom)
        Atoms in layer
    lattice : Lattice
        Unit cell lattice parameters
    layerName : str
        Unique identifier for Layer object
    '''
    
    # properties ------------------------------------------------------------------------
    atoms = property(lambda self: self._atoms,
                     lambda self, val: self.setParam(atoms=val))
    
    lattice = property(lambda self: self._lattice,
                       lambda self, val: self.setParam(lattice=val))
    
    layerName = property(lambda self: self._layerName,
                         lambda self, val: self.setParam(layerName=val))
    
    # creates instance of Layer object --------------------------------------------------
    def __init__(self, atoms, lattice, layerName):
        
        # initialize parameters
        self._layerName = None
        self._atoms = None
        self._lattice = None

        self.setParam(atoms, lattice, layerName)

        return
    
    # sets layer parameters -------------------------------------------------------------
    def setParam(self, atoms=None, lattice=None, layerName=None):
        if atoms is not None:
            self._atoms = atoms
        if lattice is not None:
            self._lattice = lattice
        if layerName is not None:
            self._layerName = layerName
        
        return
    
    # prints layer name and atoms -------------------------------------------------------
    def display(self):
        lyr = "\n".join(("Layer Name: " + self.layerName, 
                       "Atoms: " + str(self.atoms)))
        print(lyr)
    
    # prints atom labels and positions ---------------------------------------------------
    def atom_info(self):
        for i in self.atoms:
            print(i.atomLabel + ", " + str(i.xyz))


#########################################################################################
''' Child object class -- Generates translated child layer from parent layer'''
#########################################################################################

class Child (object):
    '''
    Parameters
    ----------
    layerName : str
        Unique identifier for Child object
    parent : Layer
        Generator parent layer
    transVec : nparray
        Translation vector to generate child from parent
    '''
    
    # properties ------------------------------------------------------------------------
    layerName = property(lambda self: self._layerName,
                         lambda self, val: self.setParam(layerName=val))
    
    parent = property(lambda self: self._parent,
                      lambda self, val: self.setParam(parent=val))
    
    transVec = property(lambda self: self._transVec,
                        lambda self, val: self.setParam(transVec=val))
    
    lattice = property(lambda self: self._lattice,
                       lambda self, val: self.setParam(lattice=val))
    
    # creates instance of Child object --------------------------------------------------
    def __init__(self, layerName, parent, transVec):
        
        # initialize parameters
        self._layerName = None
        self._parent = None
        self._transVec = None
        
        self._atoms = []
        self._lattice = None
        
        self.setLayerName(layerName)
        
        self.generate(parent, transVec)
        
        return
    
    # sets layer parameters -------------------------------------------------------------
    def setParam(self, layerName=None, parent=None, transVec=None):
        if layerName is not None:
            self._layerName = layerName
        if parent is not None:
            self._parent = parent
            latt = parent.lattice
            self._lattice = latt
        if transVec is not None:
            self._transVec = transVec
         
        return
    
    # creates child layer by translating parent layer -----------------------------------
    def generate(self, parent, transVec):
        alist = parent.atoms
        for i in range(len(alist)):
            parentAtom = cp.deepcopy(alist[i])
            
            splitLabel = parentAtom.atomLabel.split("_")
            newLabel = splitLabel[0] + "_" + self.layerName
            
            newPos = np.add(parentAtom.xyz, transVec)
            
            childAtom = LayerAtom(parentAtom.layerName, newLabel,
                                  parentAtom.element, newPos, 
                                  parentAtom.occupancy, parentAtom.lattice)
            
            self._atoms.append(childAtom)
            return
    
    # prints layer name and atoms -------------------------------------------------------
    def display(self):
        lyr = "\n".join(("Layer Name: " + self.layerName, 
                       "Atoms: " + str(self.atoms)))
        print(lyr)
    
    # prints atom labels and positions ---------------------------------------------------
    def atom_info(self):
        alist = self.atoms
        for i in range(len(alist)):
            print(alist[i].atomLabel + ", " + str(alist[i].xyz))
            
            
#########################################################################################
''' Additional layer functions ''' 
#########################################################################################

# imports atomic parameters from csv file -----------------------------------------------
def import_csv(path, fn):
    '''
    Parameters
    ----------
    path : str
        File directory where csv is stored
    fn : str
        File name

    Returns
    -------
    df : dataframe
        Atomic parameters
    '''
    df = pd.read_csv(path + fn + ".csv")
    return df

# Pulls atoms from imported dataframe ---------------------------------------------------
def get_csv_layers(df, lattice, layerNames):
    '''
    Parameters
    ----------
    df : dataframe
        Atomic parameters imported from csv file
    lattice : Lattice
        Unit cell lattice parameters
    layerNames : list (str)
        Layer names as documented in imported csv

    Returns
    -------
    layers : list (Layer)
    '''
    
    layers = []
    for i in layerNames:
        lyr = df[df["Layer"] == i]
        
        alist = []
        for j in range(len(lyr.index)):
            layerName = df.iloc[j]["Layer"]
            atomLabel = df.iloc[j]["Atom"]
            element = df.iloc[j]["Element"]
            x = df.iloc[j]["x"]
            y = df.iloc[j]["x"]
            z = df.iloc[j]["x"]
            occ = df.iloc[j]["Occupancy"]
            
            newAtom = LayerAtom(layerName, atomLabel, element, [x,y,z], occ, lattice)
            alist.append(newAtom)    
    
        newLayer = Layer(alist, lattice, i)
        layers.append(newLayer)
    
    return layers

#########################################################################################
''' Unitcell object class -- Stores unit cell properties '''
#########################################################################################

class Unitcell(object):
    '''
    Parameters
    ----------
    name : str
        Unique identifier for Unitcell object
    layers : list (Layer)
        List of layers of the unit cell
    lattice : Lattice
        Unit cell lattice parameters
    '''
    
    # properties ------------------------------------------------------------------------
    name = property(lambda self: self._name,
                    lambda self, val: self.setParam(name=val))
    
    layers = property(lambda self: self._layers,
                      lambda self, val: self.setParam(layers=val))
    
    lattice = property(lambda self: self._lattice,
                       lambda self, val: self.setParam(lattice=val))
    
    # creates instance of Unitcell object -----------------------------------------------
    def __init__(self, name, layers, lattice):
        
        # initialize parameters
        self._name = None
        self._layers = None
        self._lattice = None
        
        self.setParam(name, layers, lattice)
        
        return
    
    # sets cell parameters --------------------------------------------------------------
    def setParam(self, name=None, layers=None, lattice=None):
        if name is not None:
            self._name = name
        if layers is not None:
            self._layers = layers
        if lattice is not None:
            self._lattice = lattice
            
        return
    
    # prints layer names ----------------------------------------------------------------
    def layer_info(self):
        layerList = self.layers
        for i in range(len(layerList)):
            print(layerList[i].layerName)
    
    # prints atom labels and positions for each layer -----------------------------------
    def atom_info(self):
        for i in self.layers:
            for j in i.atoms:
                print(i.layerName + ", " + j.atomLabel + ", " + str(j.xyz))
            
            
#########################################################################################
''' Supercell object class -- Stores supercell properties '''
#########################################################################################

class Supercell(object):
    '''
    Parameters
    ----------
    unitcell : Unitcell
        Base unit cell
    nStacks : int
        Number of stacks / unit cells in supercell
    fltLayer : str, optional
        Name of faulted layer. The default is None.
    stackVec : nparray, optional
        Stacking vector relative to ideal position. The default is None.
    stackProb : float, optional
        Stacking probability. The default is None.
    '''
    
    # properties ------------------------------------------------------------------------
    unitcell = property(lambda self: self._unitcell)
    
    lattice = property(lambda self: self._lattice)
    
    nStacks = property(lambda self: self._nStacks, 
                       lambda self, val: self.setParam(nStacks=val))
    
    layers = property(lambda self: self._layers, 
                      lambda self, val: self.setLayers(layers=val))
    
    fltLayer = property(lambda self: self._fltLayer,
                        lambda self, val: self.setParam(fltLayer=val))
    
    stackVec = property(lambda self: self._stackVec, 
                    lambda self, val: self.setParam(stackVec=val))
    
    stackProb = property(lambda self: self._stackProb, 
                    lambda self, val: self.setParam(stackProb=val))
    
    # creates instance of Supercell object ----------------------------------------------
    def __init__(self, unitcell, nStacks, fltLayer=None, stackVec=None, stackProb=None):

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
    
    # set cell parameters ---------------------------------------------------------------
    def setParam(self, nStacks=None):
        
        if nStacks is not None:
            self._nStacks = nStacks
            
        return
    
    # define layers of supercell --------------------------------------------------------
    def setLayers(self, unitcell, fltLayer=None, stackVec=None, stackProb=None):
        
        newLayers = []
        
        for lyr in self.unitcell.layers:
            
            for n in range(self.nStacks):
                
                tag = "_n" + str(n+1)
                p = r.randint(0,100)
                
                if fltLayer is not None:
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
                    newLayerName = lyr.layerName + tag + "_fault"
                    newLyr.setParam(layerName=newLayerName, lattice=self.lattice)
                    
                    for atom in newLyr.atoms:
                        alabel = atom.atomLabel.split("_")
                        
                        newXYZ = [atom.x, atom.y, ((atom.z + n) / self.nStacks)]
                        fltXYZ = np.add(newXYZ, stackVec)
                        
                        atom.setParam(layerName=newLayerName, 
                                      atomLabel=alabel[0], xyz=fltXYZ, 
                                      lattice=self.lattice)
                if isFlt == False:
                    newLayerName = lyr.layerName + tag
                    newLyr.setParam(layerName=newLayerName, lattice=self.lattice)
                    
                    for atom in newLyr.atoms:
                        alabel = atom.atomLabel.split("_")
                        
                        newXYZ = [atom.x, atom.y, ((atom.z + n) / self.nStacks)]
                        
                        atom.setParam(layerName=newLayerName, 
                                      atomLabel=alabel[0], xyz=newXYZ, 
                                      lattice=self.lattice)
                
                newLayers.append(newLyr)
        
                
        self._layers = newLayers
        
        return
    
    # prints layer names ----------------------------------------------------------------
    def layer_info(self):
        for i in self.layers:
            print(i.layerName)
            
    def show_faults(self):
        for lyr in self.layers:
            if "fault" in lyr.layerName:
                print(lyr.layerName) 
                    
                    
#########################################################################################
''' Additional cell functions '''
#########################################################################################                  

def to_cif(cell, path, fn, inclPkl=False):
    '''
    Generate cif file from unit cell or supercell object

    Parameters
    ----------
    cell : Unitcell or Supercell
        Structure to parse to cif file format
    path : str
        File directory to save cif to
    fn : str
        File name
    title : str, optional
        Title for cif metadata. The default is None.
    inclPkl : bool, optional
        Generates pickle file of cif info when True. The default is False.
    '''
    
    import pickle as p
    
    lines = []
    
    # space group info
    lines.extend([
        "%-31s %s" % ("_symmetry_space_group_name_H-M", "'P1"),
        "%-31s %s" % ("_symmetry_Int_Tables_number", "1"),
        "%-31s %s" % ("_symmetry_cell_setting", "triclinic"),
        ""])
    
    # lattice parameters
    latt = cell.lattice
    lines.extend([
        "%-31s %.6g" % ("_cell_length_a", cell.lattice.a),
        "%-31s %.6g" % ("_cell_length_b", cell.lattice.b),
        "%-31s %.6g" % ("_cell_length_c", cell.lattice.c),
        "%-31s %.6g" % ("_cell_angle_alpha", cell.lattice.alpha),
        "%-31s %.6g" % ("_cell_angle_beta", cell.lattice.beta),
        "%-31s %.6g" % ("_cell_angle_gamma", cell.lattice.gamma),
        ""])
    
    # loop info
    lines.extend([
        "loop_",
        "  _atom_site_label",
        "  _atom_site_type_symbol",
        "  _atom_site_fract_x",
        "  _atom_site_fract_y",
        "  _atom_site_fract_z",
        "  _atom_site_U_iso_or_equiv",
        "  _atom_site_adp_type",
        "  _atom_site_occupancy" ])

    # atoms
    for lyr in cell.layers:
        for a in lyr.atoms:
            label = a.atomLabel
            elem = a.element
            x = a.x
            y = a.y
            z = a.z
            occ = a.occupancy
            aline = " %-5s %-3s %11.6f %11.6f %11.6f %11.6f %-5s %.4f" % (
                label, elem, x, y, z, 2.0, "Uiso", occ)
            lines.append(aline)

    with open(path + fn + ".cif", "w") as cif:
        for i in lines:
            cif.write(i + "\n")
    cif.close()
    
    if inclPkl == True:
        with open(path + fn + ".pkl", "wb") as pkl:
            p.dump(lines, pkl)
        pkl.close()
