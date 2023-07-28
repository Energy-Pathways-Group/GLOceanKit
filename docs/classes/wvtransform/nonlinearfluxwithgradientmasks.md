---
layout: default
title: nonlinearFluxWithGradientMasks
parent: WVTransform
grand_parent: Classes
nav_order: 129
mathjax: true
---

#  nonlinearFluxWithGradientMasks

returns the flux of each coefficient as determined by the nonlinear flux operation


---

## Declaration
```matlab
 [Fp,Fm,F0] = nonlinearFluxWithGradientMasks(ApmUMask,A0UMask,ApmUxMask,A0UxMask)
```
## Parameters
+ `ApmUMask`  mask for the wave portion of $$\vec{u}$$ in $$\vec{u} \cdot \nabla \vec{u}$$
+ `A0UMask`  mask for the geostrophic portion of $$\vec{u}$$ in $$\vec{u} \cdot \nabla \vec{u}$$
+ `ApmUxMask`  mask for the wave portion of $$\nabla \vec{u}$$ in $$\vec{u} \cdot \nabla \vec{u}$$
+ `A0UxMask`  mask for the geostrophic portion of $$\nabla \vec{u}$$ in $$\vec{u} \cdot \nabla \vec{u}$$

## Returns
+ `Fp`  flux into the Ap coefficients
+ `Fm`  flux into the Am coefficients
+ `F0`  flux into the A0 coefficients

## Discussion

  The masks are applied to the coefficients Ap,Am,A0 before computing the
  nonlinear flux, $$\vec{u} \cdot \nabla \vec{u}$$. This function offers
  more fine-grained control than -nonlinearFluxWithMask.
 
  The nonlinear flux used is the unforced, invicid equations.
 
                  
