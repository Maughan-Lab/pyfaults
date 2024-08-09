### pfInputTransMatrix.py

Last updated: 08/09/2024

---
`pfInputTransMatrix(path, prob, lyrs, numStacks, fltLyr, lyrDict, atomDict, latt)`

&nbsp;&nbsp; Reads PyFaults input file for transition matrix-type supercell construction parameters

&nbsp;&nbsp; **Parameters:**

* **path** (*str*) -- input file path
* **prob** (*array_like*) -- probabilities from supercell properties DataFrame (*pfInput*)
* **lyrs** (*array_like*) -- layer objects from unit cell properties DataFrame (*pfInput*)
* **numStacks** (*int*) -- number of stacks in supercell from supercell properties DataFrame (*pfInput*)
* **fltLyr** (*str*) -- faulted layer from supercell properties DataFrame (*pfInput*)
* **lyrDict** (*dict*) -- layer properties dictionary from unit cell properties DataFrame (*pfInput*)
* **atomDict** (*dict*) -- atom properties dictionary from unit cell properties DataFrame (*pfInput*)
* **latt** (*Lattice*) -- unit cell lattice from unit cell properties DataFrame (*pfInput*)

&nbsp;&nbsp; **Returns:**

* **cells** (*array_like*) -- supercell objects generated from transition matrix parameters

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
