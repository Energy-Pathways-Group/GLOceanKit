---
layout: default
title: masksForFlowConstituents
parent: WVTransform
grand_parent: Classes
nav_order: 125
mathjax: true
---

#  masksForFlowConstituents

Returns a sets of 'masks' indicating where different solution types live in the Ap, Am, A0 matrices.


---

## Declaration
```matlab
 [ApmMask,A0Mask] = masksForFlowConstituents(flowConstituents)
```
## Parameters
+ `flowConstituents`  `WVFlowConstituent` type

## Returns
+ `ApmMask`  mask for the Ap and Am matrices
+ `A0Mask`  mask for the A0 matrix

## Discussion

  Returns a sets of 'masks' (matrices with 1s or 0s) indicating where
  different solution types live in the Ap, Am, A0 matrices.
 
  Basic usage,
  ```matlab
  [ApmMask,A0Mask] = wvm.masksForFlowConstituents(WVFlowConstituent('internalGravityWave','inertialOscillation');
  ```
  will return a mask that contains 1 at the locations of the internal
  gravity waves and inertial oscillations in the Ap/Am matrices. All other
  entries will be zero.
 
  For example, if you define ``A = ApmMask .* Ap;`` then A will contain only the
  positive frequency internal gravity solutions and half the inertial solutions.
 
          
