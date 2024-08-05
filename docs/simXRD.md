### simXRD.py

Last updated: ADD

---
`fullSim(path, cif, wl, tt_max, savePath, pw=None, bg=None)`

&nbsp;&nbsp; Simulates power X-ray diffraction pattern, returns normalized intensity values

&nbsp;&nbsp; **Parameters:**

* **path** (*str*) -- CIF file location
* **cif** (*str*) -- CIF file name
* **wl** (*float*) -- instrument wavelength (Å)
* **tt_max** (*float*) -- maximum 2θ (°)
* **savePath** (*str*) -- location to save simulation data to
* **pw** (*float*, optional) -- peak broadening term (Å<sup>-1</sup>)
* **bg** (*float*, optional) -- average of normal background




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
