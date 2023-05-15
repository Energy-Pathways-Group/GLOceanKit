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
  + [`GLNetCDFSchemaBasisFunctionKey`](/classes/netcdffile/glnetcdfschemabasisfunctionkey.html) 
  + [`GLNetCDFSchemaDomainLengthKey`](/classes/netcdffile/glnetcdfschemadomainlengthkey.html) 
  + [`GLNetCDFSchemaDomainMinimumKey`](/classes/netcdffile/glnetcdfschemadomainminimumkey.html) 
  + [`GLNetCDFSchemaGridTypeKey`](/classes/netcdffile/glnetcdfschemagridtypekey.html) 
  + [`GLNetCDFSchemaIsComplexKey`](/classes/netcdffile/glnetcdfschemaiscomplexkey.html) attributes for variables
  + [`GLNetCDFSchemaIsCoordinateVariableKey`](/classes/netcdffile/glnetcdfschemaiscoordinatevariablekey.html) attributes for coordinate variables (dimensions)
  + [`GLNetCDFSchemaIsEvenlySampledKey`](/classes/netcdffile/glnetcdfschemaisevenlysampledkey.html) 
  + [`GLNetCDFSchemaIsFrequencyDomainKey`](/classes/netcdffile/glnetcdfschemaisfrequencydomainkey.html) 
  + [`GLNetCDFSchemaIsImaginaryPartKey`](/classes/netcdffile/glnetcdfschemaisimaginarypartkey.html) 
  + [`GLNetCDFSchemaIsPeridiocKey`](/classes/netcdffile/glnetcdfschemaisperidiockey.html) 
  + [`GLNetCDFSchemaIsRealPartKey`](/classes/netcdffile/glnetcdfschemaisrealpartkey.html) 
  + [`GLNetCDFSchemaIsRowVectorKey`](/classes/netcdffile/glnetcdfschemaisrowvectorkey.html) 
  + [`GLNetCDFSchemaMutableKey`](/classes/netcdffile/glnetcdfschemamutablekey.html) 
  + [`GLNetCDFSchemaProperNameKey`](/classes/netcdffile/glnetcdfschemapropernamekey.html) 
  + [`GLNetCDFSchemaSampleIntervalKey`](/classes/netcdffile/glnetcdfschemasampleintervalkey.html) 
  + [`GLNetCDFSchemaUniqueVariableIDKey`](/classes/netcdffile/glnetcdfschemauniquevariableidkey.html) 
  + [`GLNetCDFSchemaUnitsKey`](/classes/netcdffile/glnetcdfschemaunitskey.html) attributes for variables and dimensions
  + [`GLNetCDFSchemaVersionKey`](/classes/netcdffile/glnetcdfschemaversionkey.html) 
  + [`InitializeFromExistingFile`](/classes/netcdffile/initializefromexistingfile.html) 
  + [`addAttribute`](/classes/netcdffile/addattribute.html) 
  + [`addMutableDimension`](/classes/netcdffile/addmutabledimension.html) 
  + [`addVariable`](/classes/netcdffile/addvariable.html) 
  + [`close`](/classes/netcdffile/close.html) 
  + [`concatenateVariableAlongDimension`](/classes/netcdffile/concatenatevariablealongdimension.html) 
  + [`dimensionsForDimIDs`](/classes/netcdffile/dimensionsfordimids.html) 
  + [`initComplexVariable`](/classes/netcdffile/initcomplexvariable.html) 
  + [`initVariable`](/classes/netcdffile/initvariable.html) 
  + [`netCDFTypeForData`](/classes/netcdffile/netcdftypefordata.html) 
  + [`open`](/classes/netcdffile/open.html) 
  + [`readVariables`](/classes/netcdffile/readvariables.html) 
  + [`readVariablesAtIndexAlongDimension`](/classes/netcdffile/readvariablesatindexalongdimension.html) 
  + [`setVariable`](/classes/netcdffile/setvariable.html) 
  + [`sync`](/classes/netcdffile/sync.html) 
  + [`typeStringForTypeID`](/classes/netcdffile/typestringfortypeid.html) 
+ Initialization
  + [`NetCDFFile`](/classes/netcdffile/netcdffile.html) initialize an from existing or create new file
+ Working with dimensions
  + [`addDimension`](/classes/netcdffile/adddimension.html) Adds a both a new dimension and its associated coordinate variable to the NetCDF file.
  + [`dimensionWithName`](/classes/netcdffile/dimensionwithname.html) key-value Map to retrieve a NetCDFDimension object by name
  + [`dimensions`](/classes/netcdffile/dimensions.html) array of NetCDFDimension objects
+ Properties
  + [`attributes`](/classes/netcdffile/attributes.html) key-value Map of global attributes
  + [`complexVariableWithName`](/classes/netcdffile/complexvariablewithname.html) key-value Map to retrieve a NetCDFComplexVariable object by name
  + [`complexVariables`](/classes/netcdffile/complexvariables.html) array of NetCDFComplexVariable objects
  + [`variableWithName`](/classes/netcdffile/variablewithname.html) key-value Map to retrieve a NetCDFVariable object by name
+ Accessing file properties
  + [`ncid`](/classes/netcdffile/ncid.html) file handle
  + [`path`](/classes/netcdffile/path.html) file path the NetCDF file
+ Working with variables
  + [`variables`](/classes/netcdffile/variables.html) array of NetCDFVariable objects


---