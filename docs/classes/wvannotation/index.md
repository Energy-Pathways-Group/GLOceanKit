---
layout: default
title: WVAnnotation
parent: Classes
has_children: false
has_toc: false
---

#  WVAnnotation

annotates methods, properties, operations, and variables


---

## Declaration
```matlab
 classdef WVAnnotation < handle
```

## Overview
 
  The purpose of this class is twofold. First, it lets us add
  descriptions and detailedDescriptions for any method or property both
  inline in code, and externally in a markdown file. Second, its
  subclasses support annotating units and dimensions, which is not only
  useful for help documentation, but is also necessary for adding those
  variables to NetCDF files.
 
  This class looks for a detailedDescription in a .md file with the
  same name. It will look up to one folder deep from the call site.
 
  


## Topics
+ Initialization
  + [`WVAnnotation`](/classes/wvannotation/wvannotation.html) create a new instance of WVAnnotation
+ Properties
  + [`name`](/classes/wvannotation/name.html) name of the method, property, or variable
  + [`description`](/classes/wvannotation/description.html) short description of the method, property, or variable
  + [`detailedDescription`](/classes/wvannotation/detaileddescription.html) a detailed description of the method, property, or variable


---