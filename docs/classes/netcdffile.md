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
+ [Initialization](#initialization)
+ [Global Attributes](#global-attributes)
  + [`attributes`](#attributes)
  + [`addAttribute`](#addattributenamevalue)
+ [Coordinate Dimensions](#coordinate-dimensions)
  + [`dimensions`](#dimensions)
  + [`dimensionWithName`](#dimensionwithname)
  + [`addDimension`](#dimensionvariable--adddimensionnamevaluepropertiesdimlength)
  + [`addMutableDimension`](#dimensionvariable--addmutabledimensionnameproperties)
+ [Variables](#variables)
  + [`variables`]
  + [`variablesWithName`]
  + [`complexVariables`]
  + [`complexVariableWithName`]
  + [`initVariable`]
  + [`initComplexVariable`]
  + [`setVariable`]
  + [`addVariable`]
  + [`concatenateVariableAlongDimension`]
  + [`readVariables`]
  + [`readVariablesAtIndexAlongDimension`]

## Initialization

### `self = NetCDFFile(path,overwriteExisting)`
> Open an existing file or create a new empty file.
>
> Usage
> ```matlab
> ncfile = NetCDFFile(path)
> ```
> will load an existing file (if one exists) or create a new file (if none exists).
> ```matlab
> ncfile = NetCDFFile(path,'OVERWRITE_EXISTING')
> ```
> will delete any existing file and create a new file.

## Global Attributes

### `attributes`
> A `containers.Map` type that contains the key-value pairs of all global attributes in the NetCDF file.
>
> Usage
> ```matlab
> model = ncfile.attributes('model');
> ```

### `addAttribute(name,value)`
> Adds a global attribute to the NetCDF file.
>
> Usage
> `ncfile.addAttribute('model','WaveVortexModel');`
>
> Parameters
> - `name` (key) string with the name of the property
> - `value` value


## Coordinate Dimensions

### `dimensions`
> An array of NetCDFDimension objects for each coordinate dimension defined in the NetCDF file. The dimensions order in the array should reflect the underlying dimensionID defined in the NetCDF file.
>
> Usage
> ```matlab
> dim = ncfile.dimensions(dimID+1); % get the dimension with dimID
> ```

### `dimensionWithName`
> An container.Map of NetCDFDimension objects keyed by dimension name.
>
> Usage
> ```matlab
> xDim = ncfile.dimensionWithName('x');
> ```

### `[dimension,variable] = addDimension(name,value,properties,dimLength)`
> Adds a both a new dimension and its associated coordinate variable to the NetCDF file.
> 
> Usage
> ```matlab
> x = linspace(0,10,11);
> ncfile.addDimension('x',x,[]);
> ```
>
> Parameters
> - `name` string with the name of the dimension
> - `value` array of values along that dimension, or empty
> - `properties` containers.Map containing any key-value pairs to be associated with the dimension.
> - `dimLength` 
>
> Return Value
> - `dimension` a `NetCDFDimension` object with the newly create dimension
> - `variable` a `NetCDFVariable` object with the associated coordinate variable


---

### `[dimension,variable] = addMutableDimension(name,properties)`
> Add a new mutable dimension with coordinate variable to the NetCDF file.

---




