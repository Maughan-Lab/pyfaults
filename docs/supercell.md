### supercell.py

Stores supercell properties

Last updated: ADD

---
class **pyfaults.supercell.Supercell**(builtins.object)

&nbsp;&nbsp; `Supercell(unitcell, nStacks, fltLayer=None, stackVec=None, stackProb=None)`

&nbsp;&nbsp; Stores supercell properties

&nbsp;&nbsp; **Parameters:**

* **unitcell** (*Unitcell*) -- base unit cell
* **nStacks** (*int*) -- number of unit cell stacks in supercell
* **fltLayer** (*str*, optional) -- faulted layer name
* **stackVec** (*array_like*, optional) -- displacement vector \[x,y,z\] for faulted layer in fractional coordinates
* **stackProb** (*float*, optional) -- stacking fault probability

---
Methods:

`__init__(self, unitcell, nStacks, fltLayer=None, stackVec=None, stackProb=None)`

&nbsp;&nbsp; Initialization, defines Supercell defaults

`setParam(self, nStacks=None)`

&nbsp;&nbsp; Sets number of stacks in Supercell

`setLayers(self, unitcell, fltLayer=None, stackVec=None, stackProb=None)`

&nbsp;&nbsp; Constructs supercells layers

`show_faults(self)`

&nbsp;&nbsp; Prints names of faulted layers

`toCif(self, path)`

&nbsp;&nbsp; Generates a CIF file of Layer

&nbsp;&nbsp; **Parameters:**

* **path** (*str*) -- file path

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
