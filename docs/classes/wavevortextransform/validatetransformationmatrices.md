---
layout: default
title: ValidateTransformationMatrices
parent: WaveVortexTransform
grand_parent: Classes
nav_order: 68
---

#  ValidateTransformationMatrices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


---

## Discussion

  This is S*S^{-1} and therefore returns the values in
  wave-vortex space. So, C11 represents Ap and should be 1s
  where we expected Ap solutions to exist.
 
  Maybe check that max(abs(C12(:))) is very small.
