---
layout: default
title: attributes
parent: NetCDFFile
grand_parent: Classes
nav_order: 25
mathjax: true
---

#  attributes

key-value Map of global attributes


---

## Discussion

  A `containers.Map` type that contains the key-value pairs of all
  global attributes in the NetCDF file. This is intended to be
  *read only*. If you need to add a new attribute to file, use
  [`addAttribute`](#addattribute).
 
  Usage
  ```matlab
  model = ncfile.attributes('model');
  ```
  
