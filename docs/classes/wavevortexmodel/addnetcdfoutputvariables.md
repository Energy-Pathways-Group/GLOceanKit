---
layout: default
title: addNetCDFOutputVariables
parent: WaveVortexModel
grand_parent: Classes
---

#  addNetCDFOutputVariables

Add variables to list of variables to be written to the NetCDF variable during the model run.


---

## Declaration
```matlab
 addNetCDFOutputVariables(variables)
```
## Parameters
+ `variables`  strings of variable names.

## Discussion

       
  Pass strings of WaveVortexTransform state variables of the
  same name. This must be called before using any of the
  integrate methods.
 
  ```matlab
  model.addNetCDFOutputVariables('A0','u','v');
  ```
