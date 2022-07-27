---
layout: default
title: validateTransformationMatrices
parent: WaveVortexTransform
grand_parent: Classes
nav_order: 167
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

      
  This is S*S^{-1} and therefore returns the values in
  wave-vortex space. So, C11 represents Ap and should be 1s
  where we expected Ap solutions to exist.
 
  Maybe check that max(abs(C12(:))) is very small.
