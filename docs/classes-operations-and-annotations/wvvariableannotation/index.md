---
layout: default
title: WVVariableAnnotation
parent: Operations & annotations
has_children: false
has_toc: false
mathjax: true
---

#  WVVariableAnnotation

Describes a variable computed from the WVTransform


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVVariableAnnotation < <a href="/classes/wvannotation/" title="WVAnnotation">WVAnnotation</a></code></pre></div></div>

## Overview
  
  In addition to adding a name, description and detailed description of
  a given variable, you also specify its dimensions, units, and whether
  or note it has an imaginary part. These annotations are used for both
  online documentation and for writing to NetCDF files.
 
  Setting the two properties `isVariableWithLinearTimeStep` and
  `isVariableWithNonlinearTimeStep` are important for determining
  how the variable is cached, and when it is saved to a NetCDF file.
 
  Note that as a subclass of WVAnnotation, this class looks for
  a file (name).md in the directory where it is defined another other
  subdirectories. This file is then read-in to the detailed description
  that is used on the website.
 
  


## Topics
+ Initialization
  + [`WVVariableAnnotation`](/classes-operations-and-annotations/wvvariableannotation/wvvariableannotation.html) create a new instance of WVVariableAnnotation
+ Properties
  + [`dimensions`](/classes-operations-and-annotations/wvvariableannotation/dimensions.html) ordered cell array with the names of the dimensions
  + [`isComplex`](/classes-operations-and-annotations/wvvariableannotation/iscomplex.html) boolean indicating whether or not the variable may have an imaginary part
  + [`isVariableWithLinearTimeStep`](/classes-operations-and-annotations/wvvariableannotation/isvariablewithlineartimestep.html) boolean indicating whether the variable changes value with a linear time step
  + [`isVariableWithNonlinearTimeStep`](/classes-operations-and-annotations/wvvariableannotation/isvariablewithnonlineartimestep.html) boolean indicating whether the variable changes value with a non-linear time step
  + [`modelOp`](/classes-operations-and-annotations/wvvariableannotation/modelop.html) WVOperation responsible for computing this variable
  + [`units`](/classes-operations-and-annotations/wvvariableannotation/units.html) of the variable


---