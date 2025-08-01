---
layout: default
title: degreesOfFreedomForComplexMatrix
parent: WVGeometryDoublyPeriodic
grand_parent: Classes
nav_order: 10
mathjax: true
---

#  degreesOfFreedomForComplexMatrix

a matrix with the number of degrees-of-freedom at each entry


---

## Declaration
```matlab
 dof = WVGeometryDoublyPeriodic.degreesOfFreedomForComplexMatrix(Nx,Ny);
```
## Parameters
+ `Nx`  grid points in the x-direction
+ `Ny`  grid points in the y-direction

## Returns
+ `dof`  matrix containing dof

## Discussion

  A complex valued matrix A defined on a grid of size [Nx Ny]
  would has 2*Nx*Ny degrees-of-freedom at each grid point. In
  the Fourier domain, it also has 2*Nx*Ny degrees-of-freedom at
  each grid point.
 
          
