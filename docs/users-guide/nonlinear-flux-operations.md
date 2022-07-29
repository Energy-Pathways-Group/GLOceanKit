---
layout: default
title: Nonlinear flux operations
parent: Users guide
mathjax: true
---

#  Nonlinear flux operations

A the nonlinear flux defines how energy moves between the wave-vortex coefficients---these are the nonlinear terms in the equations of motion, transformed into wave-vortex space. The most basic implementation is a freely evolving, unforced (and undamped) flux. Customized forcing and damping are implemented by creating subclasses of WVNonlinearFluxOperation.

The unforced, undamped nonlinear terms in the equations of motion are,
$$
\begin{subequations}
\begin{align}
    \textrm{uNL}\equiv& u \partial_x u + v \partial_y u + w \partial_z u \\
    \textrm{vNL}\equiv&u \partial_x v + v \partial_y v + w \partial_z v \\
    \textrm{nNL}\equiv& u \partial_x \eta + v \partial_y \eta + w \left(\partial_z \eta +\eta \partial_z \ln N^2 \right)
\end{align}
\end{subequations}
$$

which, after a [transformation into wave-vortex space](/transformations/transformations.html), define the flux coefficients,

$$
\left[\begin{array}{c}
F_+ \\
F_- \\
F_0
\end{array} \right] \equiv 
    - \mathcal{L} \left[\begin{array}{c}
\textrm{uNL} \\
\textrm{vNL} \\
\textrm{nNL} 
\end{array} \right] \psi
$$

In this notation, the nonlinear equations of motion are

$$
\partial_t \left[\begin{array}{c} \hat{A}_+  \\  \hat{A}_-  \\\hat{A}_0 \end{array}\right] = \left[\begin{array}{c}
 F_+ \\
 F_- \\
F_0
\end{array} \right]
$$

To implement this in code, you need to subclass the `WVNonlinearFluxOperation` and override the `compute` method. The simplest implementation of this operation can be found in the `BoussinesqSpatial` subclass of the `WVNonlinearFluxOperation` and is only 4 lines of code,

```matlab
 function varargout = compute(self,wvt,varargin)
     varargout = cell(1,self.nVarOut);
     uNL = wvt.u .* wvt.diffX(wvt.u)   + wvt.v .* wvt.diffY(wvt.u)   + wvt.w .*  wvt.diffZF(wvt.u);
     vNL = wvt.u .* wvt.diffX(wvt.v)   + wvt.v .* wvt.diffY(wvt.v)   + wvt.w .*  wvt.diffZF(wvt.v);
     nNL = wvt.u .* wvt.diffX(wvt.eta) + wvt.v .* wvt.diffY(wvt.eta) + wvt.w .* (wvt.diffZG(wvt.eta) + wvt.eta .* wvt.dLnN2);

     [varargout{:}] = wvt.transformUVEtaToWaveVortex(uNL,vNL,nNL,wvt.t);
end
```

This operation can now be used by the `WVTransform` to compute nonlinear fluxes, and therefore also used by the `WVModel` to time step forward (integrate) the model. If you are just using a `WVTransform`, then calling

```matlab
wvt.nonlinearFluxOperation = WVBoussinesq()
```

will cause the 
