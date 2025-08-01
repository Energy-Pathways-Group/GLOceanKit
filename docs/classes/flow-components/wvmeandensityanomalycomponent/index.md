---
layout: default
title: WVMeanDensityAnomalyComponent
has_children: false
has_toc: false
mathjax: true
parent: Flow components
grand_parent: Class documentation
nav_order: 6
---

#  WVMeanDensityAnomalyComponent

Inertial oscillation solution group


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVInertialOscillationComponent < WVFlowComponent</code></pre></div></div>

## Overview
 
  


## Topics
+ Properties
  + [`enstrophyFactorForA0`](/classes/flow-components/wvmeandensityanomalycomponent/enstrophyfactorfora0.html) returns the qgpv multiplier for the A0 coefficient matrix.
+ Masks
  + [`maskOfPrimaryModesForCoefficientMatrix`](/classes/flow-components/wvmeandensityanomalycomponent/maskofprimarymodesforcoefficientmatrix.html) returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
+ Analytical solutions
  + [`meanDensityAnomalySolution`](/classes/flow-components/wvmeandensityanomalycomponent/meandensityanomalysolution.html) return a real-valued analytical solution of the mean density anomaly mode
  + [`solutionForModeAtIndex`](/classes/flow-components/wvmeandensityanomalycomponent/solutionformodeatindex.html) return the analytical solution at this index
+ Quadratic quantities
  + [`qgpvFactorForA0`](/classes/flow-components/wvmeandensityanomalycomponent/qgpvfactorfora0.html) returns the qgpv multiplier for the coefficient matrix.
  + [`totalEnergyFactorForCoefficientMatrix`](/classes/flow-components/wvmeandensityanomalycomponent/totalenergyfactorforcoefficientmatrix.html) returns the total energy multiplier for the coefficient matrix.
+ Other
  + [`WVMeanDensityAnomalyComponent`](/classes/flow-components/wvmeandensityanomalycomponent/wvmeandensityanomalycomponent.html) Inertial oscillation solution group
  + [`meanDensityAnomalySpatialTransformCoefficients`](/classes/flow-components/wvmeandensityanomalycomponent/meandensityanomalyspatialtransformcoefficients.html) 
  + [`meanDensityAnomalySpectralTransformCoefficients`](/classes/flow-components/wvmeandensityanomalycomponent/meandensityanomalyspectraltransformcoefficients.html) 


---