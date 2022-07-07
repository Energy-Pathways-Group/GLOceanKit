---
layout: default
title: setNetCDFOutputVariables
parent: WaveVortexModel
grand_parent: Classes
nav_order: 22
---

#  setNetCDFOutputVariables

Set list of variables to be written to the NetCDF variable during the model run.


---

## Declaration
```matlab
 setNetCDFOutputVariables(variables)
```
## Parameters
+ `variables`  strings of variable names.

## Discussion

       
  Pass strings of WaveVortexTransform state variables of the
  same name. This must be called before using any of the
  integrate methods.
 
  ```matlab
  model.setNetCDFOutputVariables('A0','u','v');
  ```
