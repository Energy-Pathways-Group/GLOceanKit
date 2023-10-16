---
layout: default
title: NetCDFFile
parent: NetCDFFile
grand_parent: Classes
nav_order: 20
mathjax: true
---

#  NetCDFFile

initialize an from existing or create new file


---

## Declaration
```matlab
 ncfile = NetCDFFile(path,options)
```
## Parameters
+ `path`  path to write file
+ `shouldOverwriteExisting`  (optional) boolean indicating whether or not to overwrite an existing file at the path. Default 0.

## Discussion

  Calling,
    ncfile = NetCDFFile(path)
  will load an existing file (if one exists) or create a new
  file (if none exists).
 
    ncfile = NetCDFFile(path,shouldOverwriteExisting=1)
  will delete any existing file and create a new file.
 
          - Returns: a new NetCDFFile instance
