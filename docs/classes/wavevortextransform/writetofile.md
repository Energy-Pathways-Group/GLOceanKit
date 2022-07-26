---
layout: default
title: writeToFile
parent: WaveVortexTransform
grand_parent: Classes
nav_order: 190
mathjax: true
---

#  writeToFile

Output the `WaveVortexTransform` to file.


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

            
