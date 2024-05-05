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
  + [`WVOrthogonalSolutionGroup`](/classes/solution-groups/wvorthogonalsolutiongroup/wvorthogonalsolutiongroup.html) create a new orthogonal solution group
+ Properties
  + [`abbreviatedName`](/classes/solution-groups/wvorthogonalsolutiongroup/abbreviatedname.html) abbreviated name
  + [`camelCaseName`](/classes/solution-groups/wvorthogonalsolutiongroup/camelcasename.html) name of the flow feature
  + [`nModes`](/classes/solution-groups/wvorthogonalsolutiongroup/nmodes.html) return the number of unique modes of this type
  + [`name`](/classes/solution-groups/wvorthogonalsolutiongroup/name.html) of the flow feature
  + [`wvt`](/classes/solution-groups/wvorthogonalsolutiongroup/wvt.html) reference to the wave vortex transform
+ Valid mode indices
+ Masks
  + [`maskOfConjugateModesForCoefficientMatrix`](/classes/solution-groups/wvorthogonalsolutiongroup/maskofconjugatemodesforcoefficientmatrix.html) returns a mask indicating where the redundant (conjugate )solutions live in the requested coefficient matrix.
  + [`maskOfModesForCoefficientMatrix`](/classes/solution-groups/wvorthogonalsolutiongroup/maskofmodesforcoefficientmatrix.html) returns a mask indicating where solutions live in the requested coefficient matrix.
  + [`maskOfPrimaryModesForCoefficientMatrix`](/classes/solution-groups/wvorthogonalsolutiongroup/maskofprimarymodesforcoefficientmatrix.html) returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
+ Analytical solutions
  + [`solutionForModeAtIndex`](/classes/solution-groups/wvorthogonalsolutiongroup/solutionformodeatindex.html) return the analytical solution for the mode at this index
+ Index Gymnastics
  + [`isValidConjugateModeNumber`](/classes/solution-groups/wvorthogonalsolutiongroup/isvalidconjugatemodenumber.html) returns a boolean indicating whether (k,l,j) is a valid mode number
  + [`isValidModeNumber`](/classes/solution-groups/wvorthogonalsolutiongroup/isvalidmodenumber.html) returns a boolean indicating whether (k,l,j) is a valid mode number
  + [`isValidPrimaryModeNumber`](/classes/solution-groups/wvorthogonalsolutiongroup/isvalidprimarymodenumber.html) returns a boolean indicating whether (k,l,j) is a valid mode number


---