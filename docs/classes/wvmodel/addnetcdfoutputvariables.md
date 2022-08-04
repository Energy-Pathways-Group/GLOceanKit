---
layout: default
title: addNetCDFOutputVariables
parent: WVModel
grand_parent: Classes
nav_order: 2
mathjax: true
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

       
  Pass strings of WVTransform state variables of the
  same name. This must be called before using any of the
  integrate methods.
 
  ```matlab
  model.addNetCDFOutputVariables('A0','u','v');
  ```
