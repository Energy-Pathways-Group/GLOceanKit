---
layout: default
title: writeToFile
parent: WVTransform
grand_parent: Classes
nav_order: 197
mathjax: true
---

#  writeToFile

Output the `WVTransform` to file.


---

## Declaration
```matlab
 ncfile = writeToFile(netcdfFile,variables,options)
```
## Parameters
+ `path`  path to write file
+ `variables`  strings of variable names.
+ `shouldOverwriteExisting`  (optional) boolean indicating whether or not to overwrite an existing file at the path. Default 0. 
+ `shouldAddDefaultVariables`  (optional) boolean indicating whether or not add default variables `A0`,`Ap`,`Am`,`t`. Default 1.

## Discussion

  Writes the WVTransform instance to file, with enough information to
  re-initialize. Pass additional variables to the variable list that
  should also be written to file.
 
  Subclasses should add any necessary properties or variables to the
  variable list before calling this superclass method.
 
            
