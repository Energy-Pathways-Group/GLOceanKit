---
layout: default
title: initComplexVariable
parent: NetCDFFile
grand_parent: Classes
nav_order: 34
mathjax: true
---

#  initComplexVariable

initialize a complex-valued variable


---

## Declaration
```matlab
 complexVariable = initComplexVariable(name,dimNames,properties,ncType)
```
## Parameters
+ `name`  name of the variable (a string)
+ `dimNames`  cell array containing the dimension names
+ `properties`  (optional) `containers.Map`
+ `properties`  ncType

## Returns
+ `complexVariable`  a NetCDFComplexVariable object

## Discussion

  NetCDF does not directly work with complex variables, so this
  method manages the hassle of working with the real and
  imaginary parts separately.
 
              
