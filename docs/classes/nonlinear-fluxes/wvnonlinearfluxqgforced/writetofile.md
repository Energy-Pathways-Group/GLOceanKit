---
layout: default
title: writeToFile
parent: WVNonlinearFluxQGForced
grand_parent: Classes
nav_order: 12
mathjax: true
---

#  writeToFile

write information about the nonlinear flux operation to file


---

## Declaration
```matlab
 writeToFile(ncfile,wvt)
```
## Parameters
+ `ncfile`  NetCDFFile instance that should be written to
+ `wvt`  the WVTransform associated with the nonlinear flux

## Discussion

  To enable model restarts, all subclass should override this
  method to write variables or parameters to file necessary to
  re-initialize directly from file.
 
        
Help for WVNonlinearFluxQGForced/writeToFile is inherited from superclass WVNonlinearFluxOperation