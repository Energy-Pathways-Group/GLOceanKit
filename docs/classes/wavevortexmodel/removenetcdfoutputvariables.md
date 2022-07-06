---
layout: default
title: removeNetCDFOutputVariables
parent: WaveVortexModel
---
#  removeNetCDFOutputVariables

Remove variables from the list of variables to be written to the NetCDF variable during the model run.


---

## Declaration
```matlab
 removeNetCDFOutputVariables(variables)
```
## Parameters
+ `variables`  strings of variable names.
## Discussion

       
  Pass strings of WaveVortexTransform state variables of the
  same name. This must be called before using any of the
  integrate methods.
 
  ```matlab
  model.removeNetCDFOutputVariables('A0','u','v');
  ```
