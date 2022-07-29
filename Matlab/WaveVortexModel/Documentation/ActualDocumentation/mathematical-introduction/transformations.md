---
layout: default
title: The wave-vortex transformation
parent: Mathematical introduction
mathjax: true
---

#  The wave-vortex transformation

The total transformations to and from wave-vortex space are defined so that

$$
\label{linear-transformations}
\left[\begin{array}{c} \hat{A}_+  \\  \hat{A}_-  \\\hat{A}_0 \end{array}\right]
 = \underbrace{T_\omega^{-1} \cdot S^{-1} \cdot D}_{\mathcal{L}} \cdot 
\left[\begin{array}{c}u \\v \\ \eta \end{array}\right]
\textrm{ and }
\left[\begin{array}{c}
u \\
v \\
\eta \\
w \\
p
\end{array} \right] = \underbrace{D^{-1} \cdot S \cdot T_\omega}_{\mathcal{L}^{-1}}
\left[\begin{array}{c} \hat{A}_+  \\  \hat{A}_-  \\\hat{A}_0 \end{array}\right]
$$

where it is understood that the physical variables are functions of $$(x,y,z,t)$$ and wave-vortex coefficients are indexed with $$jkl$$ and taken to be at some reference time $$t_0$$.

The full linear transformation *from* physical variables *to* wave-vortex space is $$\mathcal{L} = T_\omega^{-1} \cdot S^{-1} \cdot D\cdot$$, with inverse $$\mathcal{L}^{-1} = D^{-1} \cdot S \cdot T_\omega$$. The three parts to this transformation are,


+ $$D : (u,v,\eta) \to (\hat{u},\hat{v},\hat{\eta})$$ discrete transforms of the physical variables,
+ $$S^{-1} : (\hat{u},\hat{v},\hat{\eta}) \to (A_+,A_-,A_0)$$ transform to wave-vortex amplitude and,
+ $$T_\omega^{-1} : \left(A_+(t),A_-(t),A_0(t)\right) \to \left(A_+(t_0),A_-(t_0),A_0(t_0)\right)$$ an unwinding of wave phases to $$t_0$$.

## Discrete transformations

The discrete transformations $$D$$ include two Fourier transformations (one in the $$x$$ direction, the other in $$y$$) followed by the projection onto the vertical modes with the $$\mathcal{F}$$ and $$\mathcal{G}$$ matrices. The discrete transforms $$D$$ and its inverse is therefore,

$$
D = 
 \left[\begin{array}{ccc}
 \mathcal{F}  & 0 & 0 \\
 0 & \mathcal{F}  & 0 \\
 0 & 0 & \mathcal{G} 
 \end{array}\right] \cdot \mathcal{DFT}_y \cdot \mathcal{DFT}_x, \qquad
D^{-1} = 
    \mathcal{DFT}_x^{-1} \cdot \mathcal{DFT}_y^{-1} \cdot 
    \left[\begin{array}{ccccc}
    \mathcal{F}^{-1} & 0 & 0 & 0 & 0 \\
    0 & \mathcal{F}^{-1}  & 0 & 0 & 0 \\
    0 & 0 & \mathcal{G}^{-1}  & 0 & 0 \\
    0 & 0 & 0 & \mathcal{G}^{-1}  & 0 \\
    0 & 0 & 0 & 0 & \mathcal{F}^{-1}
    \end{array}\right]
$$

There are thus two basic operations,

$$
\begin{align}
    \hat{u}^{klj} =& \mathcal{F} \left[\mathcal{DFT}_y \left[\mathcal{DFT}_x \left[ u(x,y,z,t) \right] \right] \right] \\
    \hat{\eta}^{klj} =& \mathcal{G} \left[\mathcal{DFT}_y \left[\mathcal{DFT}_x \left[ \eta(x,y,z,t) \right] \right] \right]
\end{align}{}
$$

which, in code, are performed with,

```matlab
    u_hat = wvt.transformFromSpatialDomainWithF(u);
    w_hat = wvt.transformFromSpatialDomainWithG(w);
```

The two inverse transformations of $$D^{-1}$$ are,

$$
\begin{align}
    u(x,y,z,t) =&  \mathcal{DFT}_x^{-1} \left[\mathcal{DFT}_y^{-1} \left[ \mathcal{F}^{-1} \left[ \hat{u}^{klj} \right] \right] \right] \\
    \eta(x,y,z,t) =&  \mathcal{DFT}_x^{-1} \left[\mathcal{DFT}_y^{-1} \left[ \mathcal{G}^{-1} \left[ \hat{\eta}^{klj} \right] \right] \right]
\end{align}{}
$$

and implemented in code with

```matlab
    u = wvt.transformToSpatialDomainWithF(u_hat);
    w = wvt.transformToSpatialDomainWithG(w_hat);
```

## Wave-vortex sorting

The second part of the transformation is the wave-vortex (S)orting which takes $$(\hat{u},\hat{v},\hat{\eta})$$ and solves for the amplitudes of the wave and vortex solutions $$(A_+,A_-,A_0)$$. This $$3 \times 3$$ matrix operation is uniquely defined for *each* $$k,l,j$$. Components of the $$S$$ and $$S^{-1}$$ transformations would typically be referred to as $$S_{11}$$, $$S_{12}$$, etc, but instead we use the slightly more informative notation where

$$
    \begin{bmatrix}
    A_p^{klj}(t) \\
    A_m^{klj}(t) \\
    A_0^{klj}(t)
    \end{bmatrix} =
    \underbrace{\begin{bmatrix}
    A_pU^{klj} & A_pV^{klj} & A_pN^{klj} \\
    A_mU^{klj} & A_mV^{klj} & A_mN^{klj} \\
    A_0U^{klj} & A_0V^{klj} & A_0N^{klj} 
    \end{bmatrix}}_{\equiv S^{-1}}
    \begin{bmatrix}
    \hat{u}^{klj} \\
    \hat{v}^{klj} \\
    \hat{\eta}^{klj}
    \end{bmatrix}
$$

with inverse,

$$
    \begin{bmatrix}
    \hat{u}^{klj} \\
    \hat{v}^{klj} \\
    \hat{\eta}^{klj} \\
    \hat{w}^{klj}
    \end{bmatrix} =
    \underbrace{
    \begin{bmatrix}
    UA_p^{klj} & UA_m^{klj} & UA_0^{klj} \\
    VA_p^{klj} & VA_m^{klj} & VA_0^{klj} \\
    NA_p^{klj} & NA_m^{klj} & NA_0^{klj} \\
    WA_p^{klj} & WA_m^{klj} & 0
    \end{bmatrix}}_{\equiv S}
    \begin{bmatrix}
    A_p^{klj}(t) \\
    A_m^{klj}(t) \\
    A_p^{klj}(t)
    \end{bmatrix}
$$

In code for example, after computing $$(\hat{u},\hat{v},\hat{\eta})$$, $$A_p^{klj}(t)$$ is found with
```matlab
    Apt = wvt.ApU.*u_hat + wvt.ApV.*v_hat + wvt.ApN.*n_hat;
```
or, in reverse, 
```matlab
    u_hat = wvt.UAp.*wvt.Apt + wvt.UAm.*wvt.Amt + wvt.UA0.*wvt.A0t
```

To be slightly less abstract, note that in the rigid lid case, when $$k^2+l^2>0$$, $$j>0$$, $$S$$ and its inverse are just the internal gravity wave and geostrophic solutions,

$$
S = 
\left[\begin{array}{ccc} \frac{k\omega - il f_0}{\omega K} & \frac{k\omega + il f_0}{\omega K} & -i \frac{g}{f_0} l \\
\frac{l \omega + i k f_0}{\omega K } & \frac{l\omega - i k f_0}{\omega K } & i\frac{g}{f_0} k  \\
- \frac{K h}{\omega} &  \frac{K h}{\omega} & 1 \\
-iKh & -iKh & 0 \\
-\rho_0 g \frac{K h}{\omega} & \rho_0 g \frac{K h}{\omega} & \rho_0 g
\end{array}\right], \qquad
S^{-1} = \left[\begin{array}{ccc}
 \frac{k \omega + i l f_0}{2\omega K} & \frac{l \omega - i k f_0}{2\omega K} & - \frac{gK}{2\omega} \\
  \frac{k \omega - i l f_0}{2\omega K} & \frac{l \omega + i k f_0}{2\omega K} &  \frac{gK}{2\omega} \\
  i \frac{l h f_0}{\omega^2}& - i \frac{k h f_0}{\omega^2} &  \frac{ f_0^2}{ \omega^2}
\end{array}\right]
$$

## Phase winding

The final step of the transformation is to wind the phases to reference time $$t_0$$. This is done using the the phase winding operator $$T_\omega$$

$$
    T_\omega = 
    \left[\begin{array}{ccc}
    e^{i\omega t} & 0 & 0 \\
    0 & e^{-i\omega t} & 0 \\
    0 & 0 & 1
    \end{array}\right], \qquad
    T_\omega^{-1} = 
    \left[\begin{array}{ccc}
    e^{-i\omega t} & 0 & 0 \\
    0 & e^{i\omega t} & 0 \\
    0 & 0 & 1
    \end{array}\right].
$$

with inverse unwinding operator $$T_\omega^{-1}$$.
