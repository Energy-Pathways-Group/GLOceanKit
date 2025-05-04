---
layout: default
title: Introduction
parent: Mathematical introduction
mathjax: true
---

#  Introduction

The wave-vortex model is a spectral model that uses the basis of wave and vortex solutions of the rotating Boussinesq equations in arbitrary stratification to rewrite the physical variables of a fluid.

The fundamental *computation* performed is a linear transformation $$\mathcal{L}$$ (with inverse $$\mathcal{L}^{-1}$$) that takes physical variables $$(u,v,w,\eta,p)$$ that are functions of $$(x,y,z,t)$$ and projects them onto wave-vortex coefficients $$(A_+,A_-,A_0)$$ indexed with $$jkl$$. Specifically,

$$
\label{linear-transformations-abstract}
\left[\begin{array}{c} \hat{A}_+  \\  \hat{A}_-  \\\hat{A}_0 \end{array}\right]
 = \mathcal{L} \cdot
\left[\begin{array}{c}u \\v \\ \eta \end{array}\right]
\textrm{ and }
\left[\begin{array}{c}
u \\
v \\
\eta \\
w \\
p
\end{array} \right] = \mathcal{L}^{-1} \cdot
\left[\begin{array}{c} \hat{A}_+  \\  \hat{A}_-  \\\hat{A}_0 \end{array}\right].
$$

In the `WVTransform` the linear transformation and its inverse are implemented as,

```matlab
[Ap,Am,A0] = wvt.transformUVEtaToWaveVortex(U,V,N);
[U,V,W,N] = wvt.transformWaveVortexToUVWEta(Ap,Am,A0,t);
```

These transformations can be unit tested by generating a random flow state, and then confirming that they are, in fact, inverses of each other. In code, this is,

```matlab
[ApIO,AmIO,ApIGW,AmIGW,A0G,A0G0,A0rhobar] = wvt.generateRandomFlowState();
Ap = ApIO + ApIGW;
Am = AmIO + AmIGW;
A0 = A0G + A0G0 + A0rhobar;
```

followed by application of the transformations,

```matlab
[U,V,W,N] = wvt.transformWaveVortexToUVWEta(Ap,Am,A0,t).
[App,Amm,A00] = wvt.transformUVEtaToWaveVortex(U,V,N,t)
```
and then confirmation that, e.g., ``Ap`` and ``App`` are the same.

The fundamental *utility* of the wave-vortex model comes from the fact that the horizontally-averaged depth-integrated total energy is orthogonal (or diagonalized) in spectral space. That is to say that,

$$
\label{total_energy_relation}
   \frac{1}{2 L_x L_y}  \int_{-D}^0 \int_0^{L_y} \int_0^{L_x} (u^2+v^2 + w^2 + N^2 \eta^2)\, dx \, dy \, dz = \sum_{jkl} \alpha_{jkl}A_+^2 + \alpha_{jkl} A_-^2 + \beta_{jkl} A_0^2
$$

where the coefficients $$\alpha_{jkl}$$ and $$\beta_{jkl}$$ are simply precomputed constant coefficients. In code, $$\alpha_{jkl}$$ is denoted the `total energy factor` and is accessed with `wvt.Apm\_TE\_factor` and similarly $$\beta_{jkl}$$ with `wvt.A0\_TE\_factor`. To check this version of Plancherel's theorem in code, you can compute the depth integrated energy with,

```matlab
[u,v,w,eta] = wvt.variableWithName('u','v','w','eta');
E = trapz(wvt.z,mean(mean( u.^2 + v.^2 + w.^2 + shiftdim(wvt.N2,-2).*eta.*eta, 1 ),2 ) )/2;
```

and confirm that it is equal to

```matlab
App = wvt.Ap; Amm = wvt.Am; A00 = wvt.A0;
E = sum(sum(sum( wvt.Apm_TE_factor.*( App.*conj(App) + Amm.*conj(Amm) ) + wvt.A0_TE_factor.*( A00.*conj(A00) ) )));
```

The primary advantage of the wave-vortex formulation is that *every* degree-of-freedom in the wave-vortex model has a physical interpretation as an energetically independent linear solution (such as a wave or geostrophic vortex). Because each degree-of-freedom is energetically independent, this means that the nonlinear equations can be viewed as simply reshuffling energy between the different linear solutions.

As of December 2021 the model can be run in two configurations: non-hydrostatic with constant stratification and hydrostatic with arbitrary stratification. Non-hydrostatic arbitrary stratification has additional complications and is still in progress.
