### analyze.py

Last updated: 10/05/2024

---
`getNormVals(q, ints)`

&nbsp;&nbsp; Gets parameters for normalization

&nbsp;&nbsp; **Parameters:**

* **q** (*array_like*) -- *Q* values (Å<sup>-1</sup>)
* * **ints** (*array_like*) -- intensity values

&nbsp;&nbsp; **Returns:**

* **intsMax** (*float*) -- maximum intensity value
* **qAtIntsMax** (*float*) -- *Q* value corresponding to maximum intensity
* **maxIndex** (*int*) -- array index of maximum intensity value
* **intsMin** (*float*) -- minimum intensity value

---
`normalizeToExpt(exptQ, exptInts, q, ints)`

&nbsp;&nbsp; Normalize diffraction pattern to experimental data

&nbsp;&nbsp; **Parameters:**

* **exptQ** (*array_like*) -- experimental *Q* values (Å<sup>-1</sup>)
* **exptInts** (*array_like*) -- experimental intensity values
* **q** (*array_like*) -- *Q* values to be normalized (Å<sup>-1</sup>)
* **ints** (*array_like*) -- intensity values to be normalized

&nbsp;&nbsp; **Returns:**

* **normInts** (*array_like*) -- normalized intensity values

---
`normalizeToCalc(nofaultQ, nofaultInts, q, ints)`

&nbsp;&nbsp; Normalize diffraction pattern to reference calculated data

&nbsp;&nbsp; **Parameters:**

* **ofaultQ** (*array_like*) -- reference calculated *Q* values (Å<sup>-1</sup>)
* **ofaultInts** (*array_like*) -- reference calculated intensity values
* **q** (*array_like*) -- *Q* values to be normalized (Å<sup>-1</sup>)
* **ints** (*array_like*) -- intensity values to be normalized

&nbsp;&nbsp; **Returns:**

* **normInts** (*array_like*) -- normalized intensity values

---
`diffCurve(q1, q2, ints1, ints2)`

&nbsp;&nbsp; Calculates a difference curve between two sets of PXRD data

&nbsp;&nbsp; **Parameters:**

* **q1** (*array_like*) -- dataset 1 *Q* values (Å<sup>-1</sup>)
* **q2** (*array_like*) -- dataset 2 *Q* values (Å<sup>-1</sup>)
* **ints1** (*array_like*) -- dataset 1 intensity values
* **ints2** (*array_like*) -- dataset 2 intensity values

&nbsp;&nbsp; **Returns:**

* **diff_q** (*array_like*) -- *Q* values of difference curve
* **diff_ints** (*array_like*) -- intensity values of difference curve

---
`r2val(q1, q2, ints1, ints2)`

&nbsp;&nbsp; Calculates R<sup>2</sup> value between two sets of PXRD data

&nbsp;&nbsp; **Parameters:**

* **q1** (*array_like*) -- dataset 1 *Q* values (Å<sup>-1</sup>)
* **q2** (*array_like*) -- dataset 2 *Q* values (Å<sup>-1</sup>)
* **ints1** (*array_like*) -- dataset 1 intensity values
* **ints2** (*array_like*) -- dataset 2 intensity values

&nbsp;&nbsp; **Returns:**

* **r2** (*float*) -- R<sup>2</sup> value

---
`diff_r2(q1, q2, ints1, ints2)`

&nbsp;&nbsp; Calculates a difference curve and R<sup>2</sup> value between two sets of PXRD data

&nbsp;&nbsp; **Parameters:**

* **q1** (*array_like*) -- dataset 1 *Q* values (Å<sup>-1</sup>)
* **q2** (*array_like*) -- dataset 2 *Q* values (Å<sup>-1</sup>)
* **ints1** (*array_like*) -- dataset 1 intensity values
* **ints2** (*array_like*) -- dataset 2 intensity values

&nbsp;&nbsp; **Returns:**

* **r2** (*float*) -- R<sup>2</sup> value
* **diff_q** (*array_like*) -- *Q* values of difference curve
* **diff_ints** (*array_like*) -- intensity values of difference curve

---
`fitDiff(diff_ints1, diff_ints2)`

&nbsp;&nbsp; Calculates difference between two difference curves

&nbsp;&nbsp; **Parameters:**

* **diff_ints1** (*array_like*) -- dataset 1 difference in intensity values
* **diff_ints2** (*array_like*) -- dataset 2 difference in intensity values

&nbsp;&nbsp; **Returns:**

* **fitDiff** (*array_like*) -- difference of differences intensities

---
`simR2vals(simPath, exptPath, exptFN, exptWL, maxTT)`

&nbsp;&nbsp; Calculates R<sup>2</sup> values for each PXRD simulation in a directory against experimental data and generates a text file

&nbsp;&nbsp; **Parameters:**

* **simPath** (*str*) -- file path of simulations directory
* **exptPath** (*str*) -- file path of experimental data
* **exptFN** (*str*) -- experimental data file name
* **exptWL** (*str*) -- instrument wavelength (Å)
* **maxTT** (*str*) -- maximum 2θ (°)

&nbsp;&nbsp; **Returns:**

* **r2vals** (*array_like*) -- list of R<sup>2</sup> values

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
