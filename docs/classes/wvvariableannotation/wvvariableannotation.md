---
layout: default
title: WVVariableAnnotation
parent: WVVariableAnnotation
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  WVVariableAnnotation

create a new instance of WVVariableAnnotation


---

## Declaration
```matlab
 variableAnnotation = WVVariableAnnotation(name,dimensions,units,description,options)
```
## Parameters
+ `name`  name of the variable
+ `dimensions`  ordered list of the dimensions, or empty cell array
+ `units`  abbreviated SI units of the variable
+ `description`  short description of the variable
+ `isComplex`  (optional) indicates whether the variable has an imaginary part (default 0)
+ `detailedDescription`  (optional) detailed description of the variable

## Returns
+ `variableAnnotation`  a new instance of WVVariableAnnotation

## Discussion

  If a markdown file of the same name is in the same directory
  or child directory, it will be loaded as the detailed
  description upon initialization.
 
                  
