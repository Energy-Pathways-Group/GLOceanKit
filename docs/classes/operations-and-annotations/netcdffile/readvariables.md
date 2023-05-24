---
layout: default
title: readVariables
parent: NetCDFFile
grand_parent: Classes
nav_order: 41
mathjax: true
---

#  readVariables

read variables from file


---

## Declaration
```matlab
 varargout = readVariables(variableNames)
```
## Parameters
+ `variableNames`  (repeating) list of variable names

## Returns
+ `varargout`  (repeating) list of variable data

## Discussion

  Pass a list of variables to read and the data will be
  returned in the same order.
  
  ```matlab
  [x,y] = ncfile.readVariables('x','y');
  ```
 
        
