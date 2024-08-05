### lattice.py

Stores unit cell lattice parameters

Last updated: ADD

---
class **pyfaults.lattice.Lattice**(builtins.object)

&nbsp;&nbsp; `Lattice(a, b, c, alpha, beta, gamma)`

&nbsp;&nbsp; Stores unit cell lattice parameters

&nbsp;&nbsp; **Parameters:**

* **a** (*float*) -- unit cell vector *a* (Å)
* **b** (*float*) -- unit cell vector *b* (Å)
* **c** (*float*) -- unit cell vector *c* (Å)
* **alpha** (*float*) -- unit cell angle *α* (°)
* **beta** (*float*) -- unit cell angle *β* (°)
* **gamma** (*float*) -- unit cell angle *γ* (°)

---
Methods:

`__init__(self, a, b, c, alpha, beta, gamma)`

&nbsp;&nbsp; Initialization, defines Lattice defaults

`setParam(self, a=None, b=None, c=None, alpha=None, beta=None, gamma=None)`

&nbsp;&nbsp; Sets parameters of Lattice

`display(self)`

&nbsp;&nbsp; Prints lattice parameters

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
