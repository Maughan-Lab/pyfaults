PyFaults is an open-source Python library designed to model stacking fault disorder in crystalline materials and qualitatively assess the characteristic selective broadening effects in powder X-ray diffraction (PXRD).

---
Version 1.1.4

Copyright 2023 Colorado School of Mines

Authors: Sinclair R. Combs, <sinclaircombs@mines.edu>; Annalise Maughan, <amaughan@mines.edu>

GitHub: <https://github.com/Maughan-Lab/pyfaults>

Citation DOI:

Last Updated: 08/05/2024

---
**Installation:**
```
python -m pip install git+https://github.com/Maughan-Lab/pyfaults.git
```

**Requirements:** Python 2.7+/3+, NumPy, pandas, matplotlib, Dans_Diffraction

---
### Package Contents

[*layerAtom* class](layerAtom.md)

[*layer* class](layer.md)

[*lattice* class](lattice.md)

[*unitcell* class](unitcell.md)

[*supercell* class](supercell.md)

[*importCSV* method](importCSV.md)

[*genSupercells* method](genSupercells.md)

[*simXRD* method](simXRD.md)

[*pfInput* method](pfInput.md)

[*pfInputGridSearch* method](pfInputGridSearch.md)

[*pfInputTransMatrix* method](pfInputTransMatrix.md)

[*simulate* methods](simulate.md)

[*analyze* methods](analyze.md)

[*gridSearch* methods](gridSearch.md)
