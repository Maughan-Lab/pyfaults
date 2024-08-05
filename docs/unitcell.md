### unitcell.py

Stores unit cell properties

Last updated: 08/05/2024

---
class **pyfaults.unitcell.Unitcell**(builtins.object)

&nbsp;&nbsp; `Unitcell(name, layers, lattice)`

&nbsp;&nbsp; Stores unit cell properties

&nbsp;&nbsp; **Parameters:**

* **name** (*str*) -- unique identifier for unit cell
* **layers** (*array_like* of *Layer*) -- list of layers in unit cell
* **lattice** (*Lattice*) -- unit cell lattice parameters

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
