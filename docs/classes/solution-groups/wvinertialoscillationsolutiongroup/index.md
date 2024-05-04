---
layout: default
title: WVInertialOscillationSolutionGroup
has_children: false
has_toc: false
mathjax: true
parent: Solution groups
grand_parent: Class documentation
nav_order: 4
---

#  WVInertialOscillationSolutionGroup

Inertial oscillation solution group


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVInertialOscillationSolutionGroup < WVOrthogonalSolutionGroup</code></pre></div></div>

## Overview
 
  


## Topics
+ Other
  + [`WVInertialOscillationSolutionGroup`](/classes/solution-groups/wvinertialoscillationsolutiongroup/wvinertialoscillationsolutiongroup.html) Inertial oscillation solution group
  + [`contains`](/classes/solution-groups/wvinertialoscillationsolutiongroup/contains.html) 
  + [`inertialOscillationSpatialTransformCoefficients`](/classes/solution-groups/wvinertialoscillationsolutiongroup/inertialoscillationspatialtransformcoefficients.html) 
+ Analytical solutions
  + [`inertialOscillationSolution`](/classes/solution-groups/wvinertialoscillationsolutiongroup/inertialoscillationsolution.html) return a real-valued analytical solution of the internal gravity wave mode
  + [`isValidModeNumber`](/classes/solution-groups/wvinertialoscillationsolutiongroup/isvalidmodenumber.html) return a boolean indicating whether (k,l,j) is a valid mode for the given coefficientMatrix
  + [`isValidPrimaryModeNumber`](/classes/solution-groups/wvinertialoscillationsolutiongroup/isvalidprimarymodenumber.html) return a boolean indicating whether (k,l,j) is a primary mode for the given coefficientMatrix
  + [`maskForCoefficientMatrix`](/classes/solution-groups/wvinertialoscillationsolutiongroup/maskforcoefficientmatrix.html) returns a mask indicating where solutions live in the requested coefficient matrix.
  + [`maskForConjugateCoefficients`](/classes/solution-groups/wvinertialoscillationsolutiongroup/maskforconjugatecoefficients.html) returns a mask indicating where the redundant (conjugate )solutions live in the requested coefficient matrix.
  + [`maskForPrimaryCoefficients`](/classes/solution-groups/wvinertialoscillationsolutiongroup/maskforprimarycoefficients.html) returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
  + [`nUniqueSolutions`](/classes/solution-groups/wvinertialoscillationsolutiongroup/nuniquesolutions.html) return the number of unique solutions of this type
  + [`uniqueSolutionAtIndex`](/classes/solution-groups/wvinertialoscillationsolutiongroup/uniquesolutionatindex.html) return the analytical solution at this index


---