---
layout: default
title: WVGeostrophicComponent
has_children: false
has_toc: false
mathjax: true
parent: Flow components
grand_parent: Class documentation
nav_order: 3
---

#  WVGeostrophicComponent

Geostrophic solution group


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVGeostrophicComponent < WVFlowComponent</code></pre></div></div>

## Overview
  FlowConstituentGroup WVGeostrophicFlowGroup
  WVInternalGravityWaveFlowGroup
  WVRigidLidFlowGroup
  OrthogonalSolutionGroup
  


## Topics
+ Properties
  + [`enstrophyFactorForA0`](/classes/flow-components/wvgeostrophiccomponent/enstrophyfactorfora0.html) returns the qgpv multiplier for the A0 coefficient matrix.
+ Analytical solutions
  + [`geostrophicSolution`](/classes/flow-components/wvgeostrophiccomponent/geostrophicsolution.html) return a real-valued analytical solution of the geostrophic mode
  + [`maskOfPrimaryModesForCoefficientMatrix`](/classes/flow-components/wvgeostrophiccomponent/maskofprimarymodesforcoefficientmatrix.html) returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
  + [`normalizeGeostrophicModeProperties`](/classes/flow-components/wvgeostrophiccomponent/normalizegeostrophicmodeproperties.html) returns properties of a geostrophic solution relative to the primary mode number
  + [`solutionForModeAtIndex`](/classes/flow-components/wvgeostrophiccomponent/solutionformodeatindex.html) return the analytical solution at this index
+ Quadratic quantities
  + [`qgpvFactorForA0`](/classes/flow-components/wvgeostrophiccomponent/qgpvfactorfora0.html) returns the qgpv multiplier for the coefficient matrix.
  + [`randomAmplitudes`](/classes/flow-components/wvgeostrophiccomponent/randomamplitudes.html) returns random amplitude for a valid flow state
  + [`totalEnergyFactorForCoefficientMatrix`](/classes/flow-components/wvgeostrophiccomponent/totalenergyfactorforcoefficientmatrix.html) returns the total energy multiplier for the coefficient matrix.
+ Other
  + [`WVGeostrophicComponent`](/classes/flow-components/wvgeostrophiccomponent/wvgeostrophiccomponent.html) Geostrophic solution group
  + [`geostrophicSpatialTransformCoefficients`](/classes/flow-components/wvgeostrophiccomponent/geostrophicspatialtransformcoefficients.html) 
  + [`geostrophicSpectralTransformCoefficients`](/classes/flow-components/wvgeostrophiccomponent/geostrophicspectraltransformcoefficients.html) 
  + [`hkeFactorForA0`](/classes/flow-components/wvgeostrophiccomponent/hkefactorfora0.html) energy factor for most modes...
  + [`peFactorForA0`](/classes/flow-components/wvgeostrophiccomponent/pefactorfora0.html) 


---