---
layout: default
title: dimensions
parent: NetCDFFile
grand_parent: Classes
nav_order: 32
mathjax: true
---

#  dimensions

array of NetCDFDimension objects


---

## Discussion

  An array of NetCDFDimension objects for each coordinate dimension
  defined in the NetCDF file. The dimensions order in the array
  should reflect the underlying dimensionID defined in the NetCDF
  file.
 
  Usage
  ```matlab
  dim = ncfile.dimensions(dimID+1); % get the dimension with dimID
  ```
  
