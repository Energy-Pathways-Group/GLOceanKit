---
layout: default
title: setNetCDFOutputVariables
parent: WVModel
grand_parent: Classes
nav_order: 21
mathjax: true
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

       
  Pass strings of WVTransform state variables of the
  same name. This must be called before using any of the
  integrate methods.
 
  ```matlab
  model.setNetCDFOutputVariables('A0','u','v');
  ```
