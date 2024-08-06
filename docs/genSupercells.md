### genSupercells.py

Last updated: 08/05/2024

---
`genSupercells(unitcell, nStacks, fltLayer, probList, sVecList)`

&nbsp;&nbsp; Generates all supercell models as CIFs within a given parameter space, currently displacement-type construction only

&nbsp;&nbsp; **Parameters:**

* **unitcell** (*Unitcell*) -- base unit cell
* **nStacks** (*int*) -- number of unit cell stacks in supercell
* **fltLayer** (*str*) -- faulted layer name
* **probList** (*array_like*) -- list of stacking fault probabilities
* **sVecList** (*array_like*) -- list of displacement vectors \[x,y,z\] in fractional coordinates

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
