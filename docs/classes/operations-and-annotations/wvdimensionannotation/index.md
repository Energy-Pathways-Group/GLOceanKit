---
layout: default
title: WVDimensionAnnotation
has_children: false
has_toc: false
mathjax: true
parent: Operations & annotations
grand_parent: Class documentation
nav_order: 3
---

#  WVDimensionAnnotation

Describes a coordinate dimension of the WVTransform


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVDimensionAnnotation < <a href="/classes/wvannotation/" title="WVAnnotation">WVAnnotation</a></code></pre></div></div>

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
  + [`units`](/classes/wvdimensionannotation/units.html) of the dimension


---