---
layout: default
title: addDimensionAnnotations
parent: WVTransform
grand_parent: Classes
nav_order: 57
mathjax: true
---

#  addDimensionAnnotations

add one or more WVDimensions


---

## Declaration
```matlab
 addDimensionAnnotations(dimensionAnnotation)
```
## Parameters
+ `dimensionAnnotation`  one or more WVDimensionAnnotation objects

## Discussion

  Dimension annotations are used to described the coordinate dimensions
  used in the WVTransform. These annotations are used to add variables to
  NetCDF files, and also used to generate the online documentation.
 
  In general, users should not need to add dimensions. The default
  dimensions are added by the WVTransform upon initialization.
  Usage:
 
  ```matlab
  dimensions = WVDimensionAnnotation('x', 'm', 'x coordinate');
  dimensions.attributes('standard_name') = 'projection_x_coordinate';
  dimensions.attributes('axis') = 'X';
 
  wvt.addDimensionAnnotations(dimension);
  ```
      
