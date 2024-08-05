layerAtom.py

Stores atomic properties

By Sinclair R. Combs
Copyright 2023 Colorado School of Mines

Last updated: ADD

---
class **pyfaults.layerAtom.LayerAtom**(builtins.object)

&nbsp; `LayerAtom(layerName, atomLabel, element, xyz, occupancy, biso, lattice)`

&nbsp; Stores atomic properties

&nbsp; **Parameters:**

&nbsp; - **layerName** (*str*) -- name of layer that atom resides in
&nbsp; - **atomLabel** (*str*) -- unique identifier for atomic position
&nbsp; - **element** (*str*) -- elemental species and oxidation state
&nbsp; - **xyz** (*array_like*) -- atomic position in fractional coordinates
&nbsp; - **occupancy** (*float*) -- site occupancy (0 to 1)
&nbsp; - **biso** (*float*) -- isotropic atomic displacement parameter
&nbsp; - **lattice** (*Lattice*) -- unit cell lattice parameters

---
Methods:

`__init__(self, layerName, atomLabel, element, xyz, occupancy, biso, lattice)`

Initialization, defines LayerAtom defaults

`setParam(self, layerName=None, atomLabel=None, element=None, xyz=None, lattice=None, occupancy=None, biso=None)`

Sets parameters of LayerAtom
