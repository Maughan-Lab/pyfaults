### PyFaults
Version 1.1.4

Copyright 2023 Colorado School of Mines

Authors: Sinclair R. Combs, <sinclaircombs@mines.edu>; Annalise Maughan, <amaughan@mines.edu>

GitHub: <https://github.com/Maughan-Lab/pyfaults>

Citation DOI:

Last Updated: 06/30/2024

PyFaults is an open-source Python library designed to model stacking fault disorder in crystalline materials and qualitatively assess the characteristic selective broadening effects in powder X-ray diffraction (PXRD).

**Installation:**
```
python -m pip install git+https://github.com/Maughan-Lab/pyfaults.git
```

**Requirements:** Python 2.7+/3+, NumPy, pandas, matplotlib, Dans_Diffraction

### Package Contents

layerAtom
layer
lattice
unitcell
unitCellFromCSV
supercell
supercellAdjZ
TMsupercell
genSupercells

simXRD
simulate
gridSearch
analyze
