### layer.py

Stores unit cell layer properties

Last updated: ADD

---
class **pyfaults.layer.Layer**(builtins.object)

&nbsp;&nbsp; `Layer(atoms, lattice, layerName)`

&nbsp;&nbsp; Stores unit cell layer properties

&nbsp;&nbsp; **Parameters:**

* **atoms** (*array_like*) -- list of atoms in layer (LayerAtom objects)
* **lattice** (*Lattice*) -- unit cell lattice parameters
* **layerName** (*str*) -- unique identifier for layer

---
Methods:

`__init__(self, atoms, lattice, layerName)`

Initialization, defines Layer defaults

`setParam(self, atoms=None, lattice=None, layerName=None)`

Sets parameters of Layer

`display(self)`

Prints atoms in Layer

`toCif(self, path)`

Generates a CIF file of Layer

&nbsp;&nbsp; **Parameters:**

* **path** (*str*) -- file path

`genChildLayer(self, childName, transVec)`

Creates a new Layer object that is a copy of the original Layer displaced by a given translation vector

&nbsp;&nbsp; **Parameters:**

* **childName** (*str*) -- unique identifier for layer
* **transVec** (*array_like*) -- translation vector [x, y, z] relative to original/parent layer

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
