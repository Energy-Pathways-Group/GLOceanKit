---
layout: default
title: NetCDFFile
has_children: false
has_toc: false
mathjax: true
parent: Class documentation
nav_order: 6
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
  + [`NetCDFFile`](/classes/netcdffile/netcdffile.html) initialize an from existing or create new file
+ Accessing file properties
  + [`ncid`](/classes/netcdffile/ncid.html) file handle
  + [`path`](/classes/netcdffile/path.html) file path the NetCDF file
+ Working with dimensions
  + [`addDimension`](/classes/netcdffile/adddimension.html) Adds a both a new dimension and its associated coordinate variable to the NetCDF file.
  + [`dimensionWithName`](/classes/netcdffile/dimensionwithname.html) key-value Map to retrieve a NetCDFDimension object by name
  + [`dimensions`](/classes/netcdffile/dimensions.html) array of NetCDFDimension objects
  + [`dimensionsForDimIDs`](/classes/netcdffile/dimensionsfordimids.html) return the dimension IDs given the dimension names
+ Working with variables
  + [`addVariable`](/classes/netcdffile/addvariable.html) add a new (real or complex) variable to the file
  + [`complexVariableWithName`](/classes/netcdffile/complexvariablewithname.html) key-value Map to retrieve a NetCDFComplexVariable object by name
  + [`complexVariables`](/classes/netcdffile/complexvariables.html) array of NetCDFComplexVariable objects
  + [`concatenateVariableAlongDimension`](/classes/netcdffile/concatenatevariablealongdimension.html) append new data to an existing variable
  + [`initComplexVariable`](/classes/netcdffile/initcomplexvariable.html) initialize a complex-valued variable
  + [`initVariable`](/classes/netcdffile/initvariable.html) initialize a real-valued variable
  + [`readVariables`](/classes/netcdffile/readvariables.html) read variables from file
  + [`readVariablesAtIndexAlongDimension`](/classes/netcdffile/readvariablesatindexalongdimension.html) read variables from file at a particular index (e.g., time)
  + [`setVariable`](/classes/netcdffile/setvariable.html) add data for a variable with a given name
  + [`variableWithName`](/classes/netcdffile/variablewithname.html) key-value Map to retrieve a NetCDFVariable object by name
  + [`variables`](/classes/netcdffile/variables.html) array of NetCDFVariable objects
+ Working with global attributes
  + [`addAttribute`](/classes/netcdffile/addattribute.html) add a global attribute to the file
  + [`attributes`](/classes/netcdffile/attributes.html) key-value Map of global attributes
+ Schema keys
  + [`GLNetCDFSchemaUnitsKey`](/classes/netcdffile/glnetcdfschemaunitskey.html) Units of the variable or dimension
  + Dimensions
    + [`GLNetCDFSchemaBasisFunctionKey`](/classes/netcdffile/glnetcdfschemabasisfunctionkey.html) What basis function describe this dimension
    + [`GLNetCDFSchemaDomainLengthKey`](/classes/netcdffile/glnetcdfschemadomainlengthkey.html) The length of the domain
    + [`GLNetCDFSchemaDomainMinimumKey`](/classes/netcdffile/glnetcdfschemadomainminimumkey.html) The minimum value of the domain
    + [`GLNetCDFSchemaGridTypeKey`](/classes/netcdffile/glnetcdfschemagridtypekey.html) type of grid
    + [`GLNetCDFSchemaIsCoordinateVariableKey`](/classes/netcdffile/glnetcdfschemaiscoordinatevariablekey.html) A Boolean value that indicates whether the dimension is associated with a coordinate variable
    + [`GLNetCDFSchemaIsEvenlySampledKey`](/classes/netcdffile/glnetcdfschemaisevenlysampledkey.html) A Boolean value that indicates whether the dimension has even sampling
    + [`GLNetCDFSchemaIsFrequencyDomainKey`](/classes/netcdffile/glnetcdfschemaisfrequencydomainkey.html) A Boolean value that indicates whether the dimension is considered in the frequency (spectral) domain
    + [`GLNetCDFSchemaIsPeridiocKey`](/classes/netcdffile/glnetcdfschemaisperidiockey.html) A Boolean value that indicates whether the dimension is periodic
    + [`GLNetCDFSchemaMutableKey`](/classes/netcdffile/glnetcdfschemamutablekey.html) A Boolean value that indicates whether the dimension is mutable
    + [`GLNetCDFSchemaSampleIntervalKey`](/classes/netcdffile/glnetcdfschemasampleintervalkey.html) sample interval of the domain, if it is evenly sampled
  + Variables
    + [`GLNetCDFSchemaIsComplexKey`](/classes/netcdffile/glnetcdfschemaiscomplexkey.html) A Boolean value that indicates whether the variable is complex valued
    + [`GLNetCDFSchemaIsImaginaryPartKey`](/classes/netcdffile/glnetcdfschemaisimaginarypartkey.html) A Boolean value that indicates whether this is the complex part of the variable
    + [`GLNetCDFSchemaIsRealPartKey`](/classes/netcdffile/glnetcdfschemaisrealpartkey.html) A Boolean value that indicates whether this is the real part of the variable
    + [`GLNetCDFSchemaIsRowVectorKey`](/classes/netcdffile/glnetcdfschemaisrowvectorkey.html) A Boolean value that indicates whether the variable was defined as a row vector
    + [`GLNetCDFSchemaProperNameKey`](/classes/netcdffile/glnetcdfschemapropernamekey.html) Human readable name of the variable
    + [`GLNetCDFSchemaUniqueVariableIDKey`](/classes/netcdffile/glnetcdfschemauniquevariableidkey.html) A custom unique variable ID
+ Other
  + [`CreateNewFile`](/classes/netcdffile/createnewfile.html) 
  + [`GLNetCDFSchemaVersionKey`](/classes/netcdffile/glnetcdfschemaversionkey.html) - Topic: Schema keys
  + [`InitializeFromExistingFile`](/classes/netcdffile/initializefromexistingfile.html) 
  + [`addMutableDimension`](/classes/netcdffile/addmutabledimension.html) 
  + [`close`](/classes/netcdffile/close.html) - Topic: Accessing file properties
  + [`dump`](/classes/netcdffile/dump.html) 
  + [`netCDFTypeForData`](/classes/netcdffile/netcdftypefordata.html) 
  + [`open`](/classes/netcdffile/open.html) - Topic: Accessing file properties
  + [`sync`](/classes/netcdffile/sync.html) - Topic: Accessing file properties
  + [`typeStringForTypeID`](/classes/netcdffile/typestringfortypeid.html) 


---