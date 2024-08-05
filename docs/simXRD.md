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

&nbsp;&nbsp; **Returns:**

* **q** (*array_like*) -- calculated *Q* values (Å<sup>-1</sup>)
* **ints** (*array_like*) -- calculated intensity values, normalized

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
