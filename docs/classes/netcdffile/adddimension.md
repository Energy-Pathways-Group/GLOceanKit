---
layout: default
title: addDimension
parent: NetCDFFile
grand_parent: Classes
nav_order: 22
mathjax: true
---

#  addDimension

Adds a both a new dimension and its associated coordinate variable to the NetCDF file.


---

## Declaration
```matlab
 [dimension,variable] = addDimension(name,data,properties,dimLength)
```
## Parameters
+ `name`  string with the name of the dimension
+ `data`  array of values along that dimension, or empty
+ `properties`  containers.Map containing any key-value pairs to be associated with the dimension.
+ `dimLength`  (optional) length of the dimension

## Returns
+ `dimension`  a NetCDFDimension object with the newly create dimension
+ `variable`  a NetCDFVariable object with the associated coordinate variable

## Discussion

  Usage
  ```matlab
  x = linspace(0,10,11);
  ncfile.addDimension('x',x,[]);
  ```
 
                
