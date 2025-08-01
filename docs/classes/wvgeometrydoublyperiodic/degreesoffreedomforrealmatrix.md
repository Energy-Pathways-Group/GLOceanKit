---
layout: default
title: degreesOfFreedomForRealMatrix
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 11
mathjax: true
---

#  degreesOfFreedomForRealMatrix

a matrix with the number of degrees-of-freedom at each entry


---

## Declaration
```matlab
 matrix = WVGeometryDoublyPeriodic.degreesOfFreedomForFourierCoefficients(Nx,Ny,conjugateDimension);
```
## Parameters
+ `Nx`  grid points in the x-direction
+ `Ny`  grid points in the y-direction
+ `conjugateDimension`  (optional) set which dimension in the DFT grid is assumed to have the redundant conjugates (1 or 2), default is 2

## Returns
+ `dof`  matrix containing dof

## Discussion

  A real-valued matrix A defined on a grid of size [Nx Ny] has
  Nx*Ny degrees of freedom, one at each grid point. In the
  Fourier domain these degrees-of-freedom are more complicated,
  because some modes are strictly real-valued (k=l=0 and
  Nyquist), while others are complex, and there are redundant
  Hermitian conjugates.
 
            
