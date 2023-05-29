---
layout: default
title: NetCDFFile
parent: Operations & annotations
has_children: false
has_toc: false
mathjax: true
---

#  NetCDFFile

A class for reading and writing to NetCDF files


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef NetCDFFile < handle</code></pre></div></div>

## Overview
 
  NetCDF files are a standard file format for reading and writing data.
  This class is designed to simplify the task of adding new dimensions,
  variables, and attributes to a NetCDF file compared to using the
  built-in `ncread` and `ncwrite` functions.
 
  ```matlab
  ncfile = NetCDFFile('myfile.nc')
 
  % create two new dimensions and add them to the file
  x = linspace(0,10,11);
  y = linspace(-10,0,11);
  ncfile.addDimension('x',x);
  ncfile.addDimension('y',y);
 
  % Create new multi-dimensional variables, and add those to the file
  [X,Y] = ncgrid(x,y);
  ncfile.addVariable(X,{'x','y'});
  ncfile.addVariable(Y,{'x','y'});
  ```
 
                 
  


## Topics
+ Initializing
  + [`NetCDFFile`](/classes-operations-and-annotations/netcdffile/netcdffile.html) initialize an from existing or create new file
+ Accessing file properties
  + [`ncid`](/classes-operations-and-annotations/netcdffile/ncid.html) file handle
  + [`path`](/classes-operations-and-annotations/netcdffile/path.html) file path the NetCDF file
+ Working with dimensions
  + [`addDimension`](/classes-operations-and-annotations/netcdffile/adddimension.html) Adds a both a new dimension and its associated coordinate variable to the NetCDF file.
  + [`dimensionWithName`](/classes-operations-and-annotations/netcdffile/dimensionwithname.html) key-value Map to retrieve a NetCDFDimension object by name
  + [`dimensions`](/classes-operations-and-annotations/netcdffile/dimensions.html) array of NetCDFDimension objects
  + [`dimensionsForDimIDs`](/classes-operations-and-annotations/netcdffile/dimensionsfordimids.html) return the dimension IDs given the dimension names
+ Working with variables
  + [`addVariable`](/classes-operations-and-annotations/netcdffile/addvariable.html) add a new (real or complex) variable to the file
  + [`complexVariableWithName`](/classes-operations-and-annotations/netcdffile/complexvariablewithname.html) key-value Map to retrieve a NetCDFComplexVariable object by name
  + [`complexVariables`](/classes-operations-and-annotations/netcdffile/complexvariables.html) array of NetCDFComplexVariable objects
  + [`concatenateVariableAlongDimension`](/classes-operations-and-annotations/netcdffile/concatenatevariablealongdimension.html) append new data to an existing variable
  + [`initComplexVariable`](/classes-operations-and-annotations/netcdffile/initcomplexvariable.html) initialize a complex-valued variable
  + [`initVariable`](/classes-operations-and-annotations/netcdffile/initvariable.html) initialize a real-valued variable
  + [`readVariables`](/classes-operations-and-annotations/netcdffile/readvariables.html) read variables from file
  + [`readVariablesAtIndexAlongDimension`](/classes-operations-and-annotations/netcdffile/readvariablesatindexalongdimension.html) read variables from file at a particular index (e.g., time)
  + [`setVariable`](/classes-operations-and-annotations/netcdffile/setvariable.html) add data for a variable with a given name
  + [`variableWithName`](/classes-operations-and-annotations/netcdffile/variablewithname.html) key-value Map to retrieve a NetCDFVariable object by name
  + [`variables`](/classes-operations-and-annotations/netcdffile/variables.html) array of NetCDFVariable objects
+ Working with global attributes
  + [`addAttribute`](/classes-operations-and-annotations/netcdffile/addattribute.html) add a global attribute to the file
  + [`attributes`](/classes-operations-and-annotations/netcdffile/attributes.html) key-value Map of global attributes
+ Schema keys
  + [`GLNetCDFSchemaUnitsKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemaunitskey.html) Units of the variable or dimension
  + Dimensions
    + [`GLNetCDFSchemaBasisFunctionKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemabasisfunctionkey.html) What basis function describe this dimension
    + [`GLNetCDFSchemaDomainLengthKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemadomainlengthkey.html) The length of the domain
    + [`GLNetCDFSchemaDomainMinimumKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemadomainminimumkey.html) The minimum value of the domain
    + [`GLNetCDFSchemaGridTypeKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemagridtypekey.html) type of grid
    + [`GLNetCDFSchemaIsCoordinateVariableKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemaiscoordinatevariablekey.html) A Boolean value that indicates whether the dimension is associated with a coordinate variable
    + [`GLNetCDFSchemaIsEvenlySampledKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemaisevenlysampledkey.html) A Boolean value that indicates whether the dimension has even sampling
    + [`GLNetCDFSchemaIsFrequencyDomainKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemaisfrequencydomainkey.html) A Boolean value that indicates whether the dimension is considered in the frequency (spectral) domain
    + [`GLNetCDFSchemaIsPeridiocKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemaisperidiockey.html) A Boolean value that indicates whether the dimension is periodic
    + [`GLNetCDFSchemaMutableKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemamutablekey.html) A Boolean value that indicates whether the dimension is mutable
    + [`GLNetCDFSchemaSampleIntervalKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemasampleintervalkey.html) sample interval of the domain, if it is evenly sampled
  + Variables
    + [`GLNetCDFSchemaIsComplexKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemaiscomplexkey.html) A Boolean value that indicates whether the variable is complex valued
    + [`GLNetCDFSchemaIsImaginaryPartKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemaisimaginarypartkey.html) A Boolean value that indicates whether this is the complex part of the variable
    + [`GLNetCDFSchemaIsRealPartKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemaisrealpartkey.html) A Boolean value that indicates whether this is the real part of the variable
    + [`GLNetCDFSchemaIsRowVectorKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemaisrowvectorkey.html) A Boolean value that indicates whether the variable was defined as a row vector
    + [`GLNetCDFSchemaProperNameKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemapropernamekey.html) Human readable name of the variable
    + [`GLNetCDFSchemaUniqueVariableIDKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemauniquevariableidkey.html) A custom unique variable ID
+ Other
  + [`CreateNewFile`](/classes-operations-and-annotations/netcdffile/createnewfile.html) 
  + [`GLNetCDFSchemaVersionKey`](/classes-operations-and-annotations/netcdffile/glnetcdfschemaversionkey.html) - Topic: Schema keys
  + [`InitializeFromExistingFile`](/classes-operations-and-annotations/netcdffile/initializefromexistingfile.html) 
  + [`addMutableDimension`](/classes-operations-and-annotations/netcdffile/addmutabledimension.html) 
  + [`close`](/classes-operations-and-annotations/netcdffile/close.html) - Topic: Accessing file properties
  + [`netCDFTypeForData`](/classes-operations-and-annotations/netcdffile/netcdftypefordata.html) 
  + [`open`](/classes-operations-and-annotations/netcdffile/open.html) - Topic: Accessing file properties
  + [`sync`](/classes-operations-and-annotations/netcdffile/sync.html) - Topic: Accessing file properties
  + [`typeStringForTypeID`](/classes-operations-and-annotations/netcdffile/typestringfortypeid.html) 


---