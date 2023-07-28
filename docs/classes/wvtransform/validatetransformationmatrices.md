---
layout: default
title: validateTransformationMatrices
parent: WVTransform
grand_parent: Classes
nav_order: 183
mathjax: true
---

#  validateTransformationMatrices

used to confirm if $$S$$ and $$S^{-1}$$ are inverses


---

## Declaration
```matlab
 [C11,C21,C31,C12,C22,C32,C13,C23,C33] = validateTransformationMatrices(self)
```
## Discussion

      
  This is S^{-1}*S and therefore returns the values in
  wave-vortex space. So, C11 represents Ap and should be 1s
  where we expected Ap solutions to exist.
 
  Note that the current version of the code does NOT precondition the
  transformation matrices. Using the version that did (literally just scale
  by the F, G coefficients) produces much higher numerical accuracy.
