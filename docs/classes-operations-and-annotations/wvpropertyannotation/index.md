---
layout: default
title: WVPropertyAnnotation
parent: Operations & annotations
has_children: false
has_toc: false
mathjax: true
---

#  WVPropertyAnnotation

Describes a property of the WVTransform


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVPropertyAnnotation < <a href="/classes/wvannotation/" title="WVAnnotation">WVAnnotation</a></code></pre></div></div>

## Overview
 
  In addition to adding a name, description and detailed description of
  a given property, you can also specify the properties dimensions,
  its units, and whether it is a complex number or not. These
  annotations are used for both online documentation and for writing to
  NetCDF files.
 
  Note that as a subclass of WVAnnotation, this class looks for
  a file (name).md in the directory where it is defined another other
  subdirectories. This file is then read-in to the detailed description
  that is used on the website.
 
  


## Topics
+ Initialization
  + [`WVPropertyAnnotation`](/classes-operations-and-annotations/wvpropertyannotation/wvpropertyannotation.html) create a new instance of WVPropertyAnnotation
+ Properties
  + [`dimensions`](/classes-operations-and-annotations/wvpropertyannotation/dimensions.html) ordered cell array with the names of the dimensions
  + [`isComplex`](/classes-operations-and-annotations/wvpropertyannotation/iscomplex.html) boolean indicating whether or not the property may have an imaginary part
  + [`units`](/classes-operations-and-annotations/wvpropertyannotation/units.html) of the dimension


---