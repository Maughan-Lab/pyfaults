### layer.py

Stores unit cell layer properties

Last updated: 08/05/2024

---
class **pyfaults.layer.Layer**(builtins.object)

&nbsp;&nbsp; `Layer(atoms, lattice, layerName)`

&nbsp;&nbsp; Stores unit cell layer properties

&nbsp;&nbsp; **Parameters:**

* **atoms** (*array_like* of *LayerAtom*) -- list of atoms in layer
* **lattice** (*Lattice*) -- unit cell lattice parameters
* **layerName** (*str*) -- unique identifier for layer

---
Methods:

`__init__(self, atoms, lattice, layerName)`

&nbsp;&nbsp; Initialization, defines Layer defaults

`setParam(self, atoms=None, lattice=None, layerName=None)`

&nbsp;&nbsp; Sets parameters of Layer

`display(self)`

&nbsp;&nbsp; Prints atoms in Layer

`toCif(self, path)`

&nbsp;&nbsp; Generates a CIF file of Layer

&nbsp;&nbsp; **Parameters:**

* **path** (*str*) -- file path

`genChildLayer(self, childName, transVec)`

&nbsp;&nbsp; Creates a new Layer object that is a copy of the original Layer displaced by a given translation vector

&nbsp;&nbsp; **Parameters:**

* **childName** (*str*) -- unique identifier for layer
* **transVec** (*array_like*) -- translation vector [x, y, z] relative to original/parent layer

---
`getLayers(df, lattice, layerNames)`

&nbsp;&nbsp; Imports layer information from a dataframe, returns list of Layer objects

&nbsp;&nbsp; **Parameters:**

* **df** (*DataFrame*) -- atomic parameters
* **lattice** (*Lattice*) -- unit cell lattice parameters
* **layerNames** (*array_like*) -- defined layer names

&nbsp;&nbsp; **Returns:**

* **layers** (*array_like*) -- list of Layer objects

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
