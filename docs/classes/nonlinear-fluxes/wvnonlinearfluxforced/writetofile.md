---
layout: default
title: writeToFile
parent: WVNonlinearFluxForced
grand_parent: Classes
nav_order: 18
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
 
        
Help for WVNonlinearFluxForced/writeToFile is inherited from superclass WVNonlinearFluxOperation
