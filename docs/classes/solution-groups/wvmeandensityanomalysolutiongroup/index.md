---
layout: default
title: WVMeanDensityAnomalySolutionGroup
has_children: false
has_toc: false
mathjax: true
parent: Solution groups
grand_parent: Class documentation
nav_order: 5
---

#  WVMeanDensityAnomalySolutionGroup

Inertial oscillation solution group


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVInertialOscillationSolutionGroup < WVOrthogonalSolutionGroup</code></pre></div></div>

## Overview
 
  


## Topics
+ Other
  + [`WVMeanDensityAnomalySolutionGroup`](/classes/solution-groups/wvmeandensityanomalysolutiongroup/wvmeandensityanomalysolutiongroup.html) Inertial oscillation solution group
  + [`contains`](/classes/solution-groups/wvmeandensityanomalysolutiongroup/contains.html) 
  + [`meanDensityAnomalySpatialTransformCoefficients`](/classes/solution-groups/wvmeandensityanomalysolutiongroup/meandensityanomalyspatialtransformcoefficients.html) 
  + [`meanDensityAnomalySpectralTransformCoefficients`](/classes/solution-groups/wvmeandensityanomalysolutiongroup/meandensityanomalyspectraltransformcoefficients.html) 
+ Analytical solutions
  + [`isValidModeNumber`](/classes/solution-groups/wvmeandensityanomalysolutiongroup/isvalidmodenumber.html) return a boolean indicating whether (k,l,j) is a valid mode for the given coefficientMatrix
  + [`isValidPrimaryModeNumber`](/classes/solution-groups/wvmeandensityanomalysolutiongroup/isvalidprimarymodenumber.html) return a boolean indicating whether (k,l,j) is a primary mode for the given coefficientMatrix
  + [`maskForCoefficientMatrix`](/classes/solution-groups/wvmeandensityanomalysolutiongroup/maskforcoefficientmatrix.html) returns a mask indicating where solutions live in the requested coefficient matrix.
  + [`maskForConjugateCoefficients`](/classes/solution-groups/wvmeandensityanomalysolutiongroup/maskforconjugatecoefficients.html) returns a mask indicating where the redundant (conjugate )solutions live in the requested coefficient matrix.
  + [`maskForPrimaryCoefficients`](/classes/solution-groups/wvmeandensityanomalysolutiongroup/maskforprimarycoefficients.html) returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
  + [`meanDensityAnomalySolution`](/classes/solution-groups/wvmeandensityanomalysolutiongroup/meandensityanomalysolution.html) return a real-valued analytical solution of the mean density anomaly mode
  + [`nUniqueSolutions`](/classes/solution-groups/wvmeandensityanomalysolutiongroup/nuniquesolutions.html) return the number of unique solutions of this type
  + [`uniqueSolutionAtIndex`](/classes/solution-groups/wvmeandensityanomalysolutiongroup/uniquesolutionatindex.html) return the analytical solution at this index


---