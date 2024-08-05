### supercell.py

Stores supercell properties

Last updated: 08/05/2024

---
class **pyfaults.supercell.Supercell**(builtins.object)

&nbsp;&nbsp; `Supercell(unitcell, nStacks, conType, fltLayer=None, stackVec=None, stackProb=None, zAdj=None, intLayer=None)`

&nbsp;&nbsp; Stores supercell properties

&nbsp;&nbsp; **Parameters:**

* **unitcell** (*Unitcell*) -- base unit cell
* **nStacks** (*int*) -- number of unit cell stacks in supercell
* **conType** (*str*) -- type of supercell construction, can be 'Displacement' or 'Intercalation'
* **fltLayer** (*str*, optional) -- faulted layer name
* **stackVec** (*array_like*, optional) -- displacement vector \[x,y,z\] for faulted layer in fractional coordinates
* **stackProb** (*float*, optional) -- stacking fault probability
* **zAdj** (*float*, optional) -- z-direction displacement
* **intLayer** (*Layer*, optional) -- added intercalation layer in unit cell

---
Methods:

`__init__(self, unitcell, nStacks, conType, fltLayer=None, stackVec=None, stackProb=None, zAdj=None, intLayer=None)`

&nbsp;&nbsp; Initialization, defines Supercell defaults

`setParam(self, nStacks=None, conType=None)`

&nbsp;&nbsp; Sets number of stacks and construction type

`setDisLayers(self, unitcell, fltLayer=None, stackVec=None, stackProb=None)`

&nbsp;&nbsp; Constructs displacement supercell layers

`setIntLayers(self, unitcell, fltLayer=None, stackProb=None, zAdj=None, intLayer=None)`

&nbsp;&nbsp; Constructs intercalation supercell layers

`show_faults(self)`

&nbsp;&nbsp; Prints names of faulted layers

`toCif(self, path)`

&nbsp;&nbsp; Generates a CIF file of Layer

&nbsp;&nbsp; **Parameters:**

* **path** (*str*) -- file path

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
