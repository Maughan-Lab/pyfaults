### gridSearch.py

Last updated: 08/05/2024

---
`stepGridSearch(pRange, sxRange, syRange)`

&nbsp;&nbsp; Generates a set of step-wise probabilities and stacking vectors

&nbsp;&nbsp; **Parameters:**

* **pRange** (*array_like*) -- list of minimum fault probability, maximum fault probability, and step size
* **sxRange** (*array_like*) -- list of minimum stacking vector x-value, maximum stacking vector x-value, and step size
* **syRange** (*array_like*) -- list of minimum stacking vector y-value, maximum stacking vector y-value, and step size

&nbsp;&nbsp; **Returns:**

* **pList** (*array_like*) -- list of probabilities
* **sList** (*array_like*) -- list of vectors

---
`randGridSearch(pRange, sxRange, syRange, numVec)`

&nbsp;&nbsp; Generates a set of step-wise probabilities and randomized stacking vectors

&nbsp;&nbsp; **Parameters:**

* **pRange** (*array_like*) -- list of minimum fault probability, maximum fault probability, and step size
* **sxRange** (*array_like*) -- list of minimum stacking vector x-value, maximum stacking vector x-value, and step size
* **syRange** (*array_like*) -- list of minimum stacking vector y-value, maximum stacking vector y-value, and step size
* **numVec** (*int*) -- number of randomized stacking vectors to generate

&nbsp;&nbsp; **Returns:**

* **pList** (*array_like*) -- list of probabilities
* **sList** (*array_like*) -- list of vectors

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
