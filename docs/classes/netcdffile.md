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

## Adding attributes, dimensions and variables

### `addAttribute(name,value)`
> Adds a global attribute to the NetCDF file.
>
> - `name` (key) string with the name of the property
> - `value` value
>
> Usage: `ncfile.addAttribute('model','WaveVortexModel');`

---

### [dimension,variable] = addDimension(name,value,properties,dimLength)
> Adds a new dimension to the NetCDF file.
>
> - `name` string with the name of the dimension
> - `value` array of values along that dimension, or empty
> - `properties` containers.Map containing any key-value pairs to be associated with the dimension.
> - `dimLength` 
>
> Usage:
> ```matlab
> x = linspace(0,10,11);
> ncfile.addDimension('x',x,[]);
> ```

---

### `[dimension,variable] = addMutableDimension(name,properties)`

> Here is a description that follows with a lot of random text Here is a description that follows with a lot of random textHere is a description that follows with a lot of random text. And, oh, hey

---

### [dimension,variable] = addMutableDimension(name,properties)

  Here is a description that follows with a lot of random text Here is a description that follows with a lot of random textHere is a description that follows with a lot of random text. And, oh, hey



