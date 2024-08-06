### TMSupercell.py

Stores supercell properties constructed from a transition matrix

Last updated: 08/05/2024

---
class **pyfaults.TMSupercell.TMSupercell**(builtins.object)

&nbsp;&nbsp; `TMSupercell(unitcell, nStacks, fltLayer, stackProb, transMatrix)`

&nbsp;&nbsp; Stores supercell properties constructed from a transition matrix

&nbsp;&nbsp; **Parameters:**

* **unitcell** (*Unitcell*) -- base unit cell
* **nStacks** (*int*) -- number of unit cell stacks in supercell
* **fltLayer** (*str*) -- faulted layer name
* **stackProb** (*float*) -- stacking fault probability
* **transMatrix** (*DataFrame*) -- transition matrix DataFrame generated from input file

---
Methods:

`__init__(self, unitcell, nStacks, fltLayer, stackProb, transMatrix)`

&nbsp;&nbsp; Initialization, defines TMSupercell defaults

`setParam(self, nStacks=None)`

&nbsp;&nbsp; Sets number of stacks

`setLayers(self, nStacks, fltLayer, stackProb, transMatrix)`

&nbsp;&nbsp; Constructs supercell layers

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
