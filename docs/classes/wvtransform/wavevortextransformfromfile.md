---
layout: default
title: waveVortexTransformFromFile
parent: WVTransform
grand_parent: Classes
nav_order: 239
mathjax: true
---

#  waveVortexTransformFromFile

Initialize a WVTransform instance from an existing file


---

## Declaration
```matlab
 wvt = waveVortexTransformFromFile(path,options)
```
## Parameters
+ `path`  path to a NetCDF file
+ `iTime`  (optional) time index to initialize from (default 1).

## Discussion

  A WVTransform instance can be recreated from a NetCDF file and .mat
  sidecar file if the default variables were save to file. For example,
 
  ```matlab
    wvt = WVTransform.waveVortexTransformFromFile('cyprus-eddy.nc',iTime=Inf);
  ```
 
  will create a WVTransform instance populated with values from the last
  time point of the file.
 
  If you intend to read more than one time point from the save file, hold
  onto the NetCDFFile instance that is returned, and then call
  [`initFromNetCDFFile`](/classes/wvtransform/initfromnetcdffile.html). This
  avoids the relatively expensive operation recreating the WVTransform, and
  simply read the appropriate data from file.
 
  See also the users guide for [reading and writing to
  file](/users-guide/reading-and-writing-to-file.html).
 
        
