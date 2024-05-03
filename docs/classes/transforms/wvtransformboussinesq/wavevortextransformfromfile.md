---
layout: default
title: waveVortexTransformFromFile
parent: WVTransformBoussinesq
grand_parent: Classes
nav_order: 63
mathjax: true
---

#  waveVortexTransformFromFile

Initialize a WVTransformHydrostatic instance from an existing file


---

## Declaration
```matlab
 wvt = waveVortexTransformFromFile(path,options)
```
## Parameters
+ `path`  path to a NetCDF file
+ `iTime`  (optional) time index to initialize from (default 1)

## Discussion

  This static method is called by WVTransform.waveVortexTransformFromFile
  and should not need to be called directly.
 
        
