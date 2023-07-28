---
layout: default
title: generateRandomFlowState
parent: WVTransform
grand_parent: Classes
nav_order: 91
mathjax: true
---

#  generateRandomFlowState

Random flow state, separated out by solution type.


---

## Declaration
```matlab
 [ApIO,AmIO,ApIGW,AmIGW,A0G,A0G0,A0rhobar] = generateRandomFlowState()
```
## Discussion

  Generate a complete set of wave-vortex coefficients with variance at all
  physically realizable solution states.
 
  This is useful for testing that the transformation matrices are really
  complete and that the energy is diagonalizable.
  
  Adding the solution types together, gives a complete state.
  Ap = ApIO + ApIGW;
  Am = AmIO + AmIGW;
  A0 = A0G + A0G0 + A0rhobar;
 
    
