---
layout: default
title: WVDimensionAnnotation
parent: Classes
has_children: false
has_toc: false
---

#  WVDimensionAnnotation

Describes a coordinate dimension of the WVTransform


---

## Declaration
```matlab
 classdef WVDimensionAnnotation < WVAnnotation
```

## Overview
 
  In addition to adding a name, description and detailed description of
  a given dimension, you also specify its units. These annotations
  are used for both online documentation and for writing to NetCDF
  files.
 
  Note that as a subclass of WVAnnotation, this class looks for
  a file (name).md in the directory where it is defined another other
  subdirectories. This file is then read-in to the detailed description
  that is used on the website.
 
  


## Topics
+ Initialization
  + [`WVDimensionAnnotation`](/classes/wvdimensionannotation/wvdimensionannotation.html) create a new instance of WVDimensionAnnotation
+ Properties
  + [`units`](/classes/wvdimensionannotation/units.html) units of the dimension
  + [`name`](/classes/wvdimensionannotation/name.html) name of the method, property, or variable
  + [`description`](/classes/wvdimensionannotation/description.html) short description of the method, property, or variable
  + [`detailedDescription`](/classes/wvdimensionannotation/detaileddescription.html) a detailed description of the method, property, or variable


---