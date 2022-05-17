---
layout: default
title: NetCDFFile
parent: Classes
---
#  NetCDFFile

Initialize from an existing file or create new file.

## Overview

NetCDF files are a standard file format for reading and writing data. This class is designed to simply the task of adding new dimensions, variables, and attributes to a NetCDF file.


## Topics

---

### Initialization

`self = NetCDFFile(path,overwriteExisting)`

Calling,
```matlab
ncfile = NetCDFFile(path)
```
will load an existing file (if one exists) or create a new file (if none exists).
```matlab
ncfile = NetCDFFile(path,'OVERWRITE_EXISTING')
```
will delete any existing file and create a new file.

### Adding attributes, dimensions and variables

### `addAttribute(name,data)`
: Here is a description that follows
: And another line

### [dimension,variable] = addDimension(name,data,properties,dimLength)
> Here is a description that follows. Here is a description that follows. Here is a description that followsHere is a description that followsHere is a description that followsHere is a description that follows
 And another line

* `[dimension,variable] = addMutableDimension(name,properties)`

  Here is a description that follows with a lot of random text Here is a description that follows with a lot of random textHere is a description that follows with a lot of random text. And, oh, hey
  
* [dimension,variable] = addMutableDimension(name,properties)

  Here is a description that follows with a lot of random text Here is a description that follows with a lot of random textHere is a description that follows with a lot of random text. And, oh, hey

addAttribute(name,data)
: Here is a description that follows
: And another line

addAttribute(name,data)
: Here is a description that follows
: And another line


