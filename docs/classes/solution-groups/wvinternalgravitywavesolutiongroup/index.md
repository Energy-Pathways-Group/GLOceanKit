---
layout: default
title: WVInternalGravityWaveSolutionGroup
has_children: false
has_toc: false
mathjax: true
parent: Solution groups
grand_parent: Class documentation
nav_order: 3
---

#  WVInternalGravityWaveSolutionGroup

Geostrophic solution group


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVGeostrophicSolutionGroup < WVOrthogonalSolutionGroup</code></pre></div></div>

## Overview
 
  


## Topics
+ Other
  + [`WVInternalGravityWaveSolutionGroup`](/classes/solution-groups/wvinternalgravitywavesolutiongroup/wvinternalgravitywavesolutiongroup.html) Geostrophic solution group
  + [`contains`](/classes/solution-groups/wvinternalgravitywavesolutiongroup/contains.html) 
  + [`internalGravityWaveSpatialTransformCoefficients`](/classes/solution-groups/wvinternalgravitywavesolutiongroup/internalgravitywavespatialtransformcoefficients.html) 
  + [`internalGravityWaveSpectralTransformCoefficients`](/classes/solution-groups/wvinternalgravitywavesolutiongroup/internalgravitywavespectraltransformcoefficients.html) 
  + [`isValidConjugateModeNumber`](/classes/solution-groups/wvinternalgravitywavesolutiongroup/isvalidconjugatemodenumber.html) 
  + [`isValidInternalGravityWaveModeNumber`](/classes/solution-groups/wvinternalgravitywavesolutiongroup/isvalidinternalgravitywavemodenumber.html) 
+ Analytical solutions
  + [`internalGravityWaveSolution`](/classes/solution-groups/wvinternalgravitywavesolutiongroup/internalgravitywavesolution.html) return a real-valued analytical solution of the internal gravity wave mode
  + [`isValidModeNumber`](/classes/solution-groups/wvinternalgravitywavesolutiongroup/isvalidmodenumber.html) return a boolean indicating whether (k,l,j) is a valid mode for the given coefficientMatrix
  + [`isValidPrimaryModeNumber`](/classes/solution-groups/wvinternalgravitywavesolutiongroup/isvalidprimarymodenumber.html) return a boolean indicating whether (k,l,j) is a primary mode for the given coefficientMatrix
  + [`maskForCoefficientMatrix`](/classes/solution-groups/wvinternalgravitywavesolutiongroup/maskforcoefficientmatrix.html) returns a mask indicating where solutions live in the requested coefficient matrix.
  + [`maskForPrimaryCoefficients`](/classes/solution-groups/wvinternalgravitywavesolutiongroup/maskforprimarycoefficients.html) returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
  + [`nUniqueSolutions`](/classes/solution-groups/wvinternalgravitywavesolutiongroup/nuniquesolutions.html) return the number of unique solutions of this type
  + [`normalizeWaveModeProperties`](/classes/solution-groups/wvinternalgravitywavesolutiongroup/normalizewavemodeproperties.html) returns properties of a internal gravity wave solutions relative to the primary mode number
  + [`uniqueSolutionAtIndex`](/classes/solution-groups/wvinternalgravitywavesolutiongroup/uniquesolutionatindex.html) return the analytical solution at this index


---