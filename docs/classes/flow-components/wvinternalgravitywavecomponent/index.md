---
layout: default
title: WVInternalGravityWaveComponent
has_children: false
has_toc: false
mathjax: true
parent: Flow components
grand_parent: Class documentation
nav_order: 4
---

#  WVInternalGravityWaveComponent

Geostrophic solution group


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVGeostrophicComponent < WVFlowComponent</code></pre></div></div>

## Overview
 
  


## Topics
+ Analytical solutions
  + [`internalGravityWaveSolution`](/classes/flow-components/wvinternalgravitywavecomponent/internalgravitywavesolution.html) return a real-valued analytical solution of the internal gravity wave mode
  + [`maskOfPrimaryModesForCoefficientMatrix`](/classes/flow-components/wvinternalgravitywavecomponent/maskofprimarymodesforcoefficientmatrix.html) returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
  + [`normalizeWaveModeProperties`](/classes/flow-components/wvinternalgravitywavecomponent/normalizewavemodeproperties.html) returns properties of a internal gravity wave solutions relative to the primary mode number
  + [`solutionForModeAtIndex`](/classes/flow-components/wvinternalgravitywavecomponent/solutionformodeatindex.html) return the analytical solution at this index
+ Quadratic quantities
  + [`totalEnergyFactorForCoefficientMatrix`](/classes/flow-components/wvinternalgravitywavecomponent/totalenergyfactorforcoefficientmatrix.html) returns the total energy multiplier for the coefficient matrix.
+ Other
  + [`WVInternalGravityWaveComponent`](/classes/flow-components/wvinternalgravitywavecomponent/wvinternalgravitywavecomponent.html) Geostrophic solution group
  + [`internalGravityWaveSpatialTransformCoefficients`](/classes/flow-components/wvinternalgravitywavecomponent/internalgravitywavespatialtransformcoefficients.html) 
  + [`internalGravityWaveSpectralTransformCoefficients`](/classes/flow-components/wvinternalgravitywavecomponent/internalgravitywavespectraltransformcoefficients.html) 
  + [`summarizeModeAtIndex`](/classes/flow-components/wvinternalgravitywavecomponent/summarizemodeatindex.html) 


---