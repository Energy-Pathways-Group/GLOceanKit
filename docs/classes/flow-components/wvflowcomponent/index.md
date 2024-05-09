---
layout: default
title: WVFlowComponent
has_children: false
has_toc: false
mathjax: true
parent: Flow components
grand_parent: Class documentation
nav_order: 1
---

#  WVFlowComponent

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
  + [`WVFlowComponent`](/classes/flow-components/wvflowcomponent/wvflowcomponent.html) create a new orthogonal solution group
+ Properties
  + [`abbreviatedName`](/classes/flow-components/wvflowcomponent/abbreviatedname.html) abbreviated name
  + [`name`](/classes/flow-components/wvflowcomponent/name.html) of the flow feature
  + [`shortName`](/classes/flow-components/wvflowcomponent/shortname.html) name of the flow feature
  + [`wvt`](/classes/flow-components/wvflowcomponent/wvt.html) reference to the wave vortex transform
+ Masks
  + [`maskA0`](/classes/flow-components/wvflowcomponent/maska0.html) returns a mask indicating where solutions live in the A0 matrix.
  + [`maskAm`](/classes/flow-components/wvflowcomponent/maskam.html) returns a mask indicating where solutions live in the Am matrix.
  + [`maskAp`](/classes/flow-components/wvflowcomponent/maskap.html) returns a mask indicating where solutions live in the Ap matrix.


---