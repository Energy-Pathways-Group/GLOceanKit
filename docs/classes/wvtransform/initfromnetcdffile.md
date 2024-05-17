---
layout: default
title: initFromNetCDFFile
parent: WVTransform
grand_parent: Classes
nav_order: 116
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
  
  This is intended to be used in conjunction with
  [`waveVortexTransformFromFile`](/classes/wvtransform/wavevortextransformfromfile.html)
  e.g.,
 
  ```matlab
  [wvt,ncfile] = WVTransform.waveVortexTransformFromFile('cyprus-eddy.nc');
  t = ncfile.readVariables('t');
  for iTime=1:length(t)
      wvt.initFromNetCDFFile(ncfile,iTime=iTime)
      // some analysis
  end
  ```
 
  Note that this method only lightly checks that you are reading from a
  file that is compatible with this transform! So be careful.
 
  See also the users guide for [reading and writing to
  file](/users-guide/reading-and-writing-to-file.html).
  
        
