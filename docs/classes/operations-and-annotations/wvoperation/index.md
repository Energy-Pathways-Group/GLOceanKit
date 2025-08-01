---
layout: default
title: WVOperation
has_children: false
has_toc: false
mathjax: true
parent: Operations & annotations
grand_parent: Class documentation
nav_order: 1
---

#  WVOperation

Perform an operation and return a variable using a WVTransform


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVOperation < handle</code></pre></div></div>

## Overview
  
  A WVOperation allows you to add new functionality to the WVTransform by
  defining one or more variables and creating an operation for computing
  those variables.
  
  If the operation is simple and can be computed in one line, you can
  directly instantiate the WVOperation class and pass a function handle.
  For example,
  
  ```matlab
  outputVar = WVVariableAnnotation('zeta_z',{'x','y','z'},'1/s^2', 'vertical component of relative vorticity');
  f = @(wvt) wvt.diffX(wvt.v) - wvt.diffY(wvt.u);
  wvt.addOperation(WVOperation('zeta_z',outputVar,f));
  ```
  
  will enable direct calls to `wvt.zeta_z` to compute the vertical
  vorticity.
  
  Note that a `WVOperation` that computes a single variable, must have the
  same name as the variable, as specified in `WVVariableAnnotation`.
 
  More involved calculations may require subclassing WVOperation and
  overriding the `compute` method. Note that every variable returned by the
  compute operation must be described with a `WVVariableAnnotation`.
  
  


## Topics
+ Initialization
  + [`WVOperation`](/classes/operations-and-annotations/wvoperation/wvoperation.html) create a new WVOperation for computing a new variable
+ Computation
  + [`compute`](/classes/operations-and-annotations/wvoperation/compute.html) the promised variable
+ Properties
  + [`f`](/classes/operations-and-annotations/wvoperation/f.html) function handle to be called when computing the operation
  + [`nVarOut`](/classes/operations-and-annotations/wvoperation/nvarout.html) number of variables returned by the computation
  + [`name`](/classes/operations-and-annotations/wvoperation/name.html) of the operation
  + [`outputVariables`](/classes/operations-and-annotations/wvoperation/outputvariables.html) array of WVVariableAnnotations describing the outputs of the computation
+ Other
  + [`detailedDescription`](/classes/operations-and-annotations/wvoperation/detaileddescription.html) 


---