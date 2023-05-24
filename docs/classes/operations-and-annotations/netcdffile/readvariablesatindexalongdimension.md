---
layout: default
title: readVariablesAtIndexAlongDimension
parent: NetCDFFile
grand_parent: Classes
nav_order: 42
mathjax: true
---

#  readVariablesAtIndexAlongDimension

read variables from file at a particular index (e.g., time)


---

## Declaration
```matlab
 varargout = readVariables(variableNames)
```
## Parameters
+ `dimName`  name of the dimension, character string
+ `index`  index at which to read the data, positive integer
+ `variableNames`  (repeating) list of variable names

## Returns
+ `varargout`  (repeating) list of variable data

## Discussion

  Pass a list of variables to read and the data will be
  returned in the same order.
 
  ```matlab
  [u,v] = ncfile.readVariablesAtIndexAlongDimension('t',100,'u','v');
  ```
 
            
