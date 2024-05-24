---
layout: default
title: dimensionAnnotationWithName
parent: WVTransform
grand_parent: Classes
nav_order: 84
mathjax: true
---

#  dimensionAnnotationWithName

retrieve a WVDimension by name


---

## Declaration
```matlab
 dimensionAnnotationWithName(name)
```
## Parameters
+ `name`  string of dimension name

## Returns
+ `dimension`  object of WVDimensionAnnotation type

## Discussion

  Dimension annotations are used to described the coordinate dimensions
  used in the WVTransform. These annotations are used to add variables to
  NetCDF files, and also used to generate the online documentation.
 
  Usage:
 
  ```matlab
  dimension = wvt.dimensionAnnotationWithName('x');
  ```
 
        
