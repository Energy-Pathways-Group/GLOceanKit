---
layout: default
title: WVPrimaryFlowComponent
has_children: false
has_toc: false
mathjax: true
parent: Flow components
grand_parent: Class documentation
nav_order: 2
---

#  WVPrimaryFlowComponent

Orthogonal solution group


---

## Overview
 
  Each degree-of-freedom in the model is associated with an analytical
  solution to the equations of motion. This class groups together
  solutions of a particular type and provides a mapping between their
  analytical solutions and their numerical representation.
 
  Perhaps the most complicate part of the numerical implementation is
  the indexing---finding where each solution is represented
  numerically. In general, a solution will have some properties, e.g.,
    (kMode,lMode,jMode,phi,A,omegasign) 
  which will have a primary and conjugate part, each of which might be
  in two different matrices.
 
            


## Topics
+ Initialization
  + [`WVPrimaryFlowComponent`](/classes/flow-components/wvprimaryflowcomponent/wvprimaryflowcomponent.html) create a new orthogonal solution group
+ Properties
  + [`nModes`](/classes/flow-components/wvprimaryflowcomponent/nmodes.html) return the number of unique modes of this type
  + [`qgpvFactorForA0`](/classes/flow-components/wvprimaryflowcomponent/qgpvfactorfora0.html) returns the qgpv multiplier for the A0 coefficient matrix.
+ Masks
  + [`maskOfConjugateModesForCoefficientMatrix`](/classes/flow-components/wvprimaryflowcomponent/maskofconjugatemodesforcoefficientmatrix.html) returns a mask indicating where the redundant (conjugate )solutions live in the requested coefficient matrix.
  + [`maskOfModesForCoefficientMatrix`](/classes/flow-components/wvprimaryflowcomponent/maskofmodesforcoefficientmatrix.html) returns a mask indicating where solutions live in the requested coefficient matrix.
  + [`maskOfPrimaryModesForCoefficientMatrix`](/classes/flow-components/wvprimaryflowcomponent/maskofprimarymodesforcoefficientmatrix.html) returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
+ Quadratic quantities
  + [`enstrophyFactorForA0`](/classes/flow-components/wvprimaryflowcomponent/enstrophyfactorfora0.html) returns the enstrophy multiplier for the A0 coefficient matrix.
  + [`randomAmplitudes`](/classes/flow-components/wvprimaryflowcomponent/randomamplitudes.html) returns random amplitude for a valid flow state
  + [`totalEnergyFactorForCoefficientMatrix`](/classes/flow-components/wvprimaryflowcomponent/totalenergyfactorforcoefficientmatrix.html) returns the total energy multiplier for the coefficient matrix.
+ Valid mode indices
+ Analytical solutions
  + [`solutionForModeAtIndex`](/classes/flow-components/wvprimaryflowcomponent/solutionformodeatindex.html) return the analytical solution for the mode at this index
+ Index Gymnastics
  + [`isValidConjugateModeNumber`](/classes/flow-components/wvprimaryflowcomponent/isvalidconjugatemodenumber.html) returns a boolean indicating whether (k,l,j) is a valid mode number
  + [`isValidModeNumber`](/classes/flow-components/wvprimaryflowcomponent/isvalidmodenumber.html) returns a boolean indicating whether (k,l,j) is a valid mode number
  + [`isValidPrimaryModeNumber`](/classes/flow-components/wvprimaryflowcomponent/isvalidprimarymodenumber.html) returns a boolean indicating whether (k,l,j) is a valid mode number


---