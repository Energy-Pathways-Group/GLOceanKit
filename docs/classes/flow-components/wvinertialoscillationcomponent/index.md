---
layout: default
title: WVInertialOscillationComponent
has_children: false
has_toc: false
mathjax: true
parent: Flow components
grand_parent: Class documentation
nav_order: 5
---

#  WVInertialOscillationComponent

Inertial oscillation solution group


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVInertialOscillationComponent < WVFlowComponent</code></pre></div></div>

## Overview
 
  


## Topics
+ Analytical solutions
  + [`inertialOscillationSolution`](/classes/flow-components/wvinertialoscillationcomponent/inertialoscillationsolution.html) return a real-valued analytical solution of the internal gravity wave mode
  + [`maskOfConjugateModesForCoefficientMatrix`](/classes/flow-components/wvinertialoscillationcomponent/maskofconjugatemodesforcoefficientmatrix.html) returns a mask indicating where the redundant (conjugate )solutions live in the requested coefficient matrix.
  + [`maskOfPrimaryModesForCoefficientMatrix`](/classes/flow-components/wvinertialoscillationcomponent/maskofprimarymodesforcoefficientmatrix.html) returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
  + [`solutionForModeAtIndex`](/classes/flow-components/wvinertialoscillationcomponent/solutionformodeatindex.html) return the analytical solution at this index
+ Quadratic quantities
  + [`randomAmplitudes`](/classes/flow-components/wvinertialoscillationcomponent/randomamplitudes.html) returns random amplitude for a valid flow state
  + [`totalEnergyFactorForCoefficientMatrix`](/classes/flow-components/wvinertialoscillationcomponent/totalenergyfactorforcoefficientmatrix.html) returns the total energy multiplier for the coefficient matrix.
+ Other
  + [`WVInertialOscillationComponent`](/classes/flow-components/wvinertialoscillationcomponent/wvinertialoscillationcomponent.html) Inertial oscillation solution group
  + [`inertialOscillationSpatialTransformCoefficients`](/classes/flow-components/wvinertialoscillationcomponent/inertialoscillationspatialtransformcoefficients.html) 


---