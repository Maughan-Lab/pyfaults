### importCSV.py

Last updated: 08/05/2024

---
`importCSV(path, fn, lattParams, lyrNames)`

&nbsp;&nbsp; Generates a Unitcell object from CSV file with atomic parameters

&nbsp;&nbsp; See [CSVexample.csv](CSVexample.csv) for formatting requirements

&nbsp;&nbsp; **Parameters:**

* **path** (*str*) -- CSV file location
* **fn** (*str*) -- CSV file name
* **lattParams** (*array_like*) -- lattice parameters \[a, b, c, α, β, γ\]
* **lyrNames** (*array_like*) -- defined layer names

&nbsp;&nbsp; **Returns:**

* **unitcell** (*Unitcell*) -- new unit cell

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
