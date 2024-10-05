### gridSearchPlot.py

Last updated: 10/05/2024

---
`gridSearchScatterPlot(r2vals, sVec, prob, inclZ=False)`

&nbsp;&nbsp; Generates a scatter plot with grid search results, R<sup>2</sup> values are color-mapped

&nbsp;&nbsp; **Parameters:**

* **r2vals** (*array_like*) -- R<sup>2</sup> values
* **sVec** (*array_like*) -- stacking vectors as [x,y,z] arrays
* **prob** (*array_like*) -- fault probabilities
* **inclZ** (*boolean*) -- sets z-axis to stacking vector z-component if False, sets z-axis to fault probability if True

&nbsp;&nbsp; **Returns:**

* **p** (*Figure*) -- matplotlib Figure object with scatter plot

---
`gridSearchSurfacePlot(r2vals, sVec, prob, inclZ=False)`

&nbsp;&nbsp; Generates a surface plot with grid search results, R<sup>2</sup> values are color-mapped

&nbsp;&nbsp; **Parameters:**

* **r2vals** (*array_like*) -- R<sup>2</sup> values
* **sVec** (*array_like*) -- stacking vectors as [x,y,z] arrays
* **prob** (*array_like*) -- fault probabilities
* **inclZ** (*boolean*) -- sets z-axis to stacking vector z-component if False, sets z-axis to fault probability if True

&nbsp;&nbsp; **Returns:**

* **p** (*Figure*) -- matplotlib Figure object with surface plot

---
By Sinclair R. Combs

Copyright 2023 Colorado School of Mines
