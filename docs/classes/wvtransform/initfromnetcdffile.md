---
layout: default
title: initFromNetCDFFile
parent: WVTransform
grand_parent: Classes
nav_order: 101
mathjax: true
---

#  initFromNetCDFFile

initialize the flow from a NetCDF file


---

## Declaration
```matlab
 initWithFile(self,path,options)
```
## Parameters
+ `ncfile`  a NetCDF file object
+ `iTime`  (optional) time index to initialize from (default 1)

## Discussion

  Clears variables Ap,Am,A0 and then sets them the values found in the file
  at the requested time.
  
  Note that this method only lightly checks that you are reading from a
  file that is compatible with this transform! So be careful
  
        
