---
layout: default
title: concatenateVariableAlongDimension
parent: NetCDFFile
grand_parent: Classes
nav_order: 30
mathjax: true
---

#  concatenateVariableAlongDimension

append new data to an existing variable


---

## Declaration
```matlab
 concatenateVariableAlongDimension(name,data,dimName,index)
```
## Parameters
+ `name`  name of the variable (a string)
+ `data`  variable data
+ `dimName`  the variable dimension along which to concatenate
+ `properties`  index at which to write data

## Discussion

  concatenates data along a variable dimension (such as a time
  dimension).
 
  ```matlab
  ncfile.concatenateVariableAlongDimension('u',u,'t',outputIndex);
  ```
 
            
