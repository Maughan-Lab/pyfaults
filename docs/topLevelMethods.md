### Top-Level Methods in __init__.py

Last updated: 08/13/2024

---
`toCif(cell, path, fn)`

&nbsp;&nbsp; Export unit cell or supercell as a CIF

&nbsp;&nbsp; **Parameters:**

* **cell** (*Unitcell* or *Supercell*) -- unit cell or supercell object
* **path** (*str*) -- CIF file location
* **fn** (*str*) -- CIF file name

---
`tt_to_q(twotheta, wavelength)`

&nbsp;&nbsp; Converts 2θ (°) values to *Q* (Å<sup>-1</sup>) values

&nbsp;&nbsp; **Parameters:**

* **twotheta** (*array_like*) -- 2θ (°) values
* **wavelength** (*float*) -- instrument wavelength (Å)

&nbsp;&nbsp; **Returns:**

* **Q** (*array_like*) -- *Q* (Å<sup>-1</sup>) values

---
`norm(ints)`

&nbsp;&nbsp; Normalizes intensity values

&nbsp;&nbsp; **Parameters:**

* **ints** (*array_like*) -- intensity values

&nbsp;&nbsp; **Returns:**

* **norm_ints** (*array_like*) -- normalized intensity values

---
`importCSV(path, fn)`

&nbsp;&nbsp; Imports a CSV file with atomic parameters

&nbsp;&nbsp; **Parameters:**

* **path** (*str*) -- CSV file path
* **fn** (*str*) -- CSV file name

&nbsp;&nbsp; **Returns:**

* **df** (*DataFrame*) -- DataFrame containing atomic parameter information

---
`importFile(path, fn, ext=None, norm=True)`

&nbsp;&nbsp; Imports X-ray diffraction data files

&nbsp;&nbsp; **Parameters:**

* **path** (*str*) -- file path
* **fn** (*str*) -- file name
* **ext** (*str*, optional) -- file extension, defaults to '.txt'
* **norm** (*bool*, optional) -- set to True to normalize intensity data, defaults to False

&nbsp;&nbsp; **Returns:**

* **q** (*array_like*) -- *Q* (Å<sup>-1</sup>) values
* **ints** (*array_like*) -- intensity values

---
`importExpt(path, fn, wl, maxTT, ext=None)`

&nbsp;&nbsp; Imports experimental X-ray diffraction data

&nbsp;&nbsp; **Parameters:**

* **path** (*str*) -- file path
* **fn** (*str*) -- file name
* **wl** (*float*) -- instrument wavelength (Å)
* **maxTT** (*float*) -- maximum 2θ (°) value
* **ext** (*str*, optional) -- file extension, defaults to '.txt'

&nbsp;&nbsp; **Returns:**

* **exptQ** (*array_like*) -- *Q* (Å<sup>-1</sup>) values
* **exptInts** (*array_like*) -- normalized intensity values

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
