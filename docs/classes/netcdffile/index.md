---
layout: default
title: NetCDFFile
parent: Classes
has_children: false
has_toc: false
mathjax: true
---

#  NetCDFFile

A class for reading and writing to NetCDF files


---



## Topics
+ Other
  + [`CreateNewFile`](/classes/netcdffile/createnewfile.html) 
  + [`GLNetCDFSchemaBasisFunctionKey`](/classes/netcdffile/glnetcdfschemabasisfunctionkey.html) - Topic: Schema keys — Dimensions
  + [`GLNetCDFSchemaDomainLengthKey`](/classes/netcdffile/glnetcdfschemadomainlengthkey.html) - Topic: Schema keys — Dimensions
  + [`GLNetCDFSchemaDomainMinimumKey`](/classes/netcdffile/glnetcdfschemadomainminimumkey.html) - Topic: Schema keys — Dimensions
  + [`GLNetCDFSchemaGridTypeKey`](/classes/netcdffile/glnetcdfschemagridtypekey.html) - Topic: Schema keys — Dimensions
  + [`GLNetCDFSchemaIsComplexKey`](/classes/netcdffile/glnetcdfschemaiscomplexkey.html) - Topic: Schema keys — Variables
  + [`GLNetCDFSchemaIsCoordinateVariableKey`](/classes/netcdffile/glnetcdfschemaiscoordinatevariablekey.html) - Topic: Schema keys — Dimensions
  + [`GLNetCDFSchemaIsEvenlySampledKey`](/classes/netcdffile/glnetcdfschemaisevenlysampledkey.html) - Topic: Schema keys — Dimensions
  + [`GLNetCDFSchemaIsFrequencyDomainKey`](/classes/netcdffile/glnetcdfschemaisfrequencydomainkey.html) - Topic: Schema keys — Dimensions
  + [`GLNetCDFSchemaIsImaginaryPartKey`](/classes/netcdffile/glnetcdfschemaisimaginarypartkey.html) - Topic: Schema keys — Variables
  + [`GLNetCDFSchemaIsPeridiocKey`](/classes/netcdffile/glnetcdfschemaisperidiockey.html) - Topic: Schema keys — Dimensions
  + [`GLNetCDFSchemaIsRealPartKey`](/classes/netcdffile/glnetcdfschemaisrealpartkey.html) - Topic: Schema keys — Variables
  + [`GLNetCDFSchemaIsRowVectorKey`](/classes/netcdffile/glnetcdfschemaisrowvectorkey.html) - Topic: Schema keys — Variables
  + [`GLNetCDFSchemaMutableKey`](/classes/netcdffile/glnetcdfschemamutablekey.html) - Topic: Schema keys — Dimensions
  + [`GLNetCDFSchemaProperNameKey`](/classes/netcdffile/glnetcdfschemapropernamekey.html) - Topic: Schema keys — Variables
  + [`GLNetCDFSchemaSampleIntervalKey`](/classes/netcdffile/glnetcdfschemasampleintervalkey.html) - Topic: Schema keys — Dimensions
  + [`GLNetCDFSchemaUniqueVariableIDKey`](/classes/netcdffile/glnetcdfschemauniquevariableidkey.html) - Topic: Schema keys — Variables
  + [`GLNetCDFSchemaUnitsKey`](/classes/netcdffile/glnetcdfschemaunitskey.html) - Topic: Schema keys
  + [`GLNetCDFSchemaVersionKey`](/classes/netcdffile/glnetcdfschemaversionkey.html) - Topic: Schema keys
  + [`InitializeFromExistingFile`](/classes/netcdffile/initializefromexistingfile.html) 
  + [`addAttribute`](/classes/netcdffile/addattribute.html) - Topic: Working with global attributes
  + [`addMutableDimension`](/classes/netcdffile/addmutabledimension.html) 
  + [`close`](/classes/netcdffile/close.html) - Topic: Accessing file properties
  + [`netCDFTypeForData`](/classes/netcdffile/netcdftypefordata.html) 
  + [`open`](/classes/netcdffile/open.html) - Topic: Accessing file properties
  + [`sync`](/classes/netcdffile/sync.html) - Topic: Accessing file properties
  + [`typeStringForTypeID`](/classes/netcdffile/typestringfortypeid.html) 
+ Initialization
  + [`NetCDFFile`](/classes/netcdffile/netcdffile.html) initialize an from existing or create new file
+ Working with dimensions
  + [`addDimension`](/classes/netcdffile/adddimension.html) Adds a both a new dimension and its associated coordinate variable to the NetCDF file.
  + [`dimensionWithName`](/classes/netcdffile/dimensionwithname.html) key-value Map to retrieve a NetCDFDimension object by name
  + [`dimensions`](/classes/netcdffile/dimensions.html) array of NetCDFDimension objects
  + [`dimensionsForDimIDs`](/classes/netcdffile/dimensionsfordimids.html) return the dimension IDs given the dimension names
+ Working with variables
  + [`addVariable`](/classes/netcdffile/addvariable.html) add a new variable to the file
  + [`complexVariableWithName`](/classes/netcdffile/complexvariablewithname.html) key-value Map to retrieve a NetCDFComplexVariable object by name
  + [`complexVariables`](/classes/netcdffile/complexvariables.html) array of NetCDFComplexVariable objects
  + [`concatenateVariableAlongDimension`](/classes/netcdffile/concatenatevariablealongdimension.html) append new data to an existing variable
  + [`initComplexVariable`](/classes/netcdffile/initcomplexvariable.html) initialize a complex-valued variable
  + [`initVariable`](/classes/netcdffile/initvariable.html) initialize a real-valued variable
  + [`readVariables`](/classes/netcdffile/readvariables.html) read data from variables
  + [`readVariablesAtIndexAlongDimension`](/classes/netcdffile/readvariablesatindexalongdimension.html) read data from variables from a particular index
  + [`setVariable`](/classes/netcdffile/setvariable.html) add data for a variable with a given name
  + [`variableWithName`](/classes/netcdffile/variablewithname.html) key-value Map to retrieve a NetCDFVariable object by name
  + [`variables`](/classes/netcdffile/variables.html) array of NetCDFVariable objects
+ Working with global attributes
  + [`attributes`](/classes/netcdffile/attributes.html) key-value Map of global attributes
+ Accessing file properties
  + [`ncid`](/classes/netcdffile/ncid.html) file handle
  + [`path`](/classes/netcdffile/path.html) file path the NetCDF file


---