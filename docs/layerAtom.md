### layerAtom.py

Stores atomic properties

Last updated: 08/05/2024

---
class **pyfaults.layerAtom.LayerAtom**(builtins.object)

&nbsp;&nbsp; `LayerAtom(layerName, atomLabel, element, xyz, occupancy, biso, lattice)`

&nbsp;&nbsp; Stores atomic properties

&nbsp;&nbsp; **Parameters:**

* **layerName** (*str*) -- name of layer that atom resides in
* **atomLabel** (*str*) -- unique identifier for atomic position
* **element** (*str*) -- elemental species and oxidation state
* **xyz** (*array_like*) -- atomic position in fractional coordinates
* **occupancy** (*float*) -- site occupancy (0 to 1)
* **biso** (*float*) -- isotropic atomic displacement parameter
* **lattice** (*Lattice*) -- unit cell lattice parameters

---
Methods:

`__init__(self, layerName, atomLabel, element, xyz, occupancy, biso, lattice)`

&nbsp;&nbsp; Initialization, defines LayerAtom defaults

`setParam(self, layerName=None, atomLabel=None, element=None, xyz=None, lattice=None, occupancy=None, biso=None)`

&nbsp;&nbsp; Sets parameters of LayerAtom

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
