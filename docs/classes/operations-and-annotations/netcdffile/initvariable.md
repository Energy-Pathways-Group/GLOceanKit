---
layout: default
title: initVariable
parent: NetCDFFile
grand_parent: Classes
nav_order: 36
mathjax: true
---

#  initVariable

initialize a real-valued variable


---

## Declaration
```matlab
 variable = initVariable(name,dimNames,properties,ncType)
```
## Parameters
+ `name`  name of the variable (a string)
+ `dimNames`  cell array containing the dimension names
+ `properties`  (optional) `containers.Map`
+ `properties`  ncType

## Returns
+ `variable`  a NetCDFVariable object

## Discussion

  The basic work flow is that you need to first,
  - initVariable followed by,
  - setVariable
  Or you can just call
  - addVariable
  which will both initialize and then set. The reason to
 
  ```matlab
  ncfile.initVariable('fluid-tracer', {'x','y','t'},containers.Map({'isTracer'},{'1'}),'NC_DOUBLE');
  ncfile.setVariable('fluid-tracer',myVariable);
  ```
 
              
