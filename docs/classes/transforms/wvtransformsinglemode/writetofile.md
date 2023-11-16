---
layout: default
title: writeToFile
parent: WVTransformSingleMode
grand_parent: Classes
nav_order: 26
mathjax: true
---

#  writeToFile

Output the `WVTransformSingleMode` instance to file.


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

            
