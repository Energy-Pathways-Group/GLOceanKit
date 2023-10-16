---
layout: default
title: addVariable
parent: NetCDFFile
grand_parent: Classes
nav_order: 24
mathjax: true
---

#  addVariable

add a new (real or complex) variable to the file


---

## Declaration
```matlab
 variable = addVariable(name,data,dimNames,properties,ncType)
```
## Parameters
+ `name`  name of the variable (a string)
+ `data`  variable data
+ `dimNames`  cell array containing the dimension names
+ `properties`  (optional) `containers.Map`
+ `properties`  ncType

## Returns
+ `variable`  a NetCDFVariable object

## Discussion

  Immediately initializes and writes the given variable data to
  file.
 
  ```matlab
  ncfile.addVariable('fluid-tracer',myVariableData, {'x','y','t'},containers.Map({'isTracer'},{'1'}),'NC_DOUBLE');
  ```
 
                
