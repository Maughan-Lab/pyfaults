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

`__init__(self, name, layers, lattice)`

&nbsp;&nbsp; Initialization, defines Unitcell defaults

`setParam(self, name=None, layers=None, lattice=None)`

&nbsp;&nbsp; Sets parameters of Unitcell

`layer_info(self)`

&nbsp;&nbsp; Prints layer names in unit cell

`toCif(self, path)`

&nbsp;&nbsp; Generates a CIF file of Layer

&nbsp;&nbsp; **Parameters:**

* **path** (*str*) -- file path

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
