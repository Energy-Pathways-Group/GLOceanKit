---
layout: default
title: WVOrthogonalSolutionGroup
parent: WVOrthogonalSolutionGroup
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  WVOrthogonalSolutionGroup

Orthogonal solution group


---

## Declaration
```matlab
 classdef WVOrthogonalSolutionGroup
```
## Discussion

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
 
  
