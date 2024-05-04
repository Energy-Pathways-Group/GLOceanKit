---
layout: default
title: WVGeostrophicSolutionGroup
has_children: false
has_toc: false
mathjax: true
parent: Solution groups
grand_parent: Class documentation
nav_order: 2
---

#  WVGeostrophicSolutionGroup

Geostrophic solution group


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVGeostrophicSolutionGroup < WVOrthogonalSolutionGroup</code></pre></div></div>

## Overview
 
  


## Topics
+ Other
  + [`WVGeostrophicSolutionGroup`](/classes/solution-groups/wvgeostrophicsolutiongroup/wvgeostrophicsolutiongroup.html) Geostrophic solution group
  + [`contains`](/classes/solution-groups/wvgeostrophicsolutiongroup/contains.html) 
  + [`geostrophicSpatialTransformCoefficients`](/classes/solution-groups/wvgeostrophicsolutiongroup/geostrophicspatialtransformcoefficients.html) 
  + [`geostrophicSpectralTransformCoefficients`](/classes/solution-groups/wvgeostrophicsolutiongroup/geostrophicspectraltransformcoefficients.html) 
  + [`isValidConjugateModeNumber`](/classes/solution-groups/wvgeostrophicsolutiongroup/isvalidconjugatemodenumber.html) Geostrophic modes are valid at all Fourier modes, except k=l=0.
  + [`isValidGeostrophicModeNumber`](/classes/solution-groups/wvgeostrophicsolutiongroup/isvalidgeostrophicmodenumber.html) 
  + [`isValidModeNumber`](/classes/solution-groups/wvgeostrophicsolutiongroup/isvalidmodenumber.html) Geostrophic modes are valid at all Fourier modes, except k=l=0.
  + [`isValidPrimaryModeNumber`](/classes/solution-groups/wvgeostrophicsolutiongroup/isvalidprimarymodenumber.html) Geostrophic modes are valid at all Fourier modes, except k=l=0.
+ Analytical solutions
  + [`geostrophicSolution`](/classes/solution-groups/wvgeostrophicsolutiongroup/geostrophicsolution.html) return a real-valued analytical solution of the geostrophic mode
  + [`maskForCoefficientMatrix`](/classes/solution-groups/wvgeostrophicsolutiongroup/maskforcoefficientmatrix.html) returns a mask indicating where solutions live in the requested coefficient matrix.
  + [`maskForPrimaryCoefficients`](/classes/solution-groups/wvgeostrophicsolutiongroup/maskforprimarycoefficients.html) returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
  + [`nUniqueSolutions`](/classes/solution-groups/wvgeostrophicsolutiongroup/nuniquesolutions.html) return the number of unique solutions of this type
  + [`normalizeGeostrophicModeProperties`](/classes/solution-groups/wvgeostrophicsolutiongroup/normalizegeostrophicmodeproperties.html) returns properties of a geostrophic solution relative to the primary mode number
  + [`uniqueSolutionAtIndex`](/classes/solution-groups/wvgeostrophicsolutiongroup/uniquesolutionatindex.html) return the analytical solution at this index


---