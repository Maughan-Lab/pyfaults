### analyze.py

Last updated: 08/09/2024

---
`pfInput(path)`

&nbsp;&nbsp; Reads PyFaults input file

&nbsp;&nbsp; **Parameters:**

* **path** (*str*) -- input file path

&nbsp;&nbsp; **Returns:**

* **unitcell** (*str*) -- unit cell
* **ucDF** (*DataFrame*) -- DataFrame with unit cell properties (structure name, simulation type, lattice, number of layers, layer names, layer dictionary, atom dictionary, layer objects)
* **gsDF** (*DataFrame*) -- DataFrame with grid search properties (type of grid search, probability range, x-component range, y-component range, number of vectors)
* **scDF** (*DataFrame*) -- DataFrame with supercell properties (faulted layer, probabilities, stacking vectors, number of stacks, z-component adjustment)
* **simDF** (*DataFrame*) -- DataFrame with simulation properties (wavelength, maximum two theta, peak broadening term)

 ---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
