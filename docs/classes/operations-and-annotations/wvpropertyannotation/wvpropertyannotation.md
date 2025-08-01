---
layout: default
title: WVPropertyAnnotation
parent: WVPropertyAnnotation
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  WVPropertyAnnotation

create a new instance of WVPropertyAnnotation


---

## Declaration
```matlab
 propAnnotation = WVPropertyAnnotation(name,dimensions,units,description,options)
```
## Parameters
+ `name`  name of the property
+ `dimensions`  ordered list of the dimensions, or empty cell array
+ `units`  abbreviated SI units of the property
+ `description`  short description of the property
+ `isComplex`  (optional) indicates whether the property has an imaginary part (default 0)
+ `detailedDescription`  (optional) detailed description of the property

## Returns
+ `propAnnotation`  a new instance of WVPropertyAnnotation

## Discussion

  If a markdown file of the same name is in the same directory
  or child directory, it will be loaded as the detailed
  description upon initialization.
 
                  
