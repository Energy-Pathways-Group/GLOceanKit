---
layout: default
title: WVOrthogonalSolutionGroup
has_children: false
has_toc: false
mathjax: true
parent: Solution groups
grand_parent: Class documentation
nav_order: 1
---

#  WVOrthogonalSolutionGroup

Orthogonal solution group


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef WVOrthogonalSolutionGroup</code></pre></div></div>

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
+ Other
  + [`WVOrthogonalSolutionGroup`](/classes/solution-groups/wvorthogonalsolutiongroup/wvorthogonalsolutiongroup.html) Orthogonal solution group
  + [`wvt`](/classes/solution-groups/wvorthogonalsolutiongroup/wvt.html) 
+ Properties
  + [`abbreviatedName`](/classes/solution-groups/wvorthogonalsolutiongroup/abbreviatedname.html) abbreviated name
  + [`camelCaseName`](/classes/solution-groups/wvorthogonalsolutiongroup/camelcasename.html) name of the flow feature
  + [`name`](/classes/solution-groups/wvorthogonalsolutiongroup/name.html) of the flow feature
+ Analytical solutions
  + [`isValidModeNumber`](/classes/solution-groups/wvorthogonalsolutiongroup/isvalidmodenumber.html) return a boolean indicating whether (k,l,j) is a valid mode for the given coefficientMatrix
  + [`isValidPrimaryModeNumber`](/classes/solution-groups/wvorthogonalsolutiongroup/isvalidprimarymodenumber.html) return a boolean indicating whether (k,l,j) is a primary mode for the given coefficientMatrix
  + [`maskForCoefficientMatrix`](/classes/solution-groups/wvorthogonalsolutiongroup/maskforcoefficientmatrix.html) returns a mask indicating where solutions live in the requested coefficient matrix.
  + [`maskForConjugateCoefficients`](/classes/solution-groups/wvorthogonalsolutiongroup/maskforconjugatecoefficients.html) returns a mask indicating where the redundant (conjugate )solutions live in the requested coefficient matrix.
  + [`maskForPrimaryCoefficients`](/classes/solution-groups/wvorthogonalsolutiongroup/maskforprimarycoefficients.html) returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
  + [`nUniqueSolutions`](/classes/solution-groups/wvorthogonalsolutiongroup/nuniquesolutions.html) return the number of unique solutions of this type
  + [`uniqueSolutionAtIndex`](/classes/solution-groups/wvorthogonalsolutiongroup/uniquesolutionatindex.html) return the analytical solution at this index


---