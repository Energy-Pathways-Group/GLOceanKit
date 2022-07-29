---
layout: default
title: Orthogonality
parent: Mathematical introduction
mathjax: true
---

#  Orthogonality in the Wave-Vortex model

The wave-vortex model is a spectral model, but unlike typical spectral models where physical variables (u,v,rho) are projected onto sines or cosines, physical variables are projected onto the linear solutions of the system, namely, the geostrophic and inertia-gravity wave solutions. A key feature that makes this model useful, is that the geostrophic and inertia-gravity wave solutions are *orthogonal* in horizontally-averaged depth-integrated total energy. As a consequence, the spectrum of this model physically meaningful.

The purpose of this document is to clarify 1) *why* orthogonality is important, and 2) *how* it is used in the model.

## Sine and cosine series

Let's first review orthogonality of sine and cosine series for both continuous and discrete system. As a simple example, let's assume we have some function in the spatial domain,

$$
f(x) = a_0 + a_1 \cos(\pi x/L) + a_2 \cos(2 \pi x/L)
$$

which has been chosen to be easily projected onto a cosine series with wavenumbers $$m_n = n \pi/L$$ where $$n=0..\infty$$. The cosine series has *orthogonality* condition,

$$
\frac{1}{L} \int_0^L \cos( n \pi x/L) \cos( m\pi x/L ) \, dx= \delta_{mn}
$$

which immediately leads to a projection operation,

$$
a_n = \frac{1}{L} \int_0^L f(x) \cos( n \pi/L) \, dx.
$$

Using the projection operation, the function $$f(x)$$ in the spatial domain can therefore also be represented as a cosine series, with coefficients $$(a_0,a_1,a_2,0,..)$$.

### Plancherel's theorem

One of the most remarkable features of the cosine series is Plancherel's theorem,

$$
\frac{1}{L} \int_0^L f^2(x) \, dx = \sum_{n=0} a_n^2.
$$

 Plancherel's theorem tells us that the *square of the sum, is the same as the sum of the squares*.

To see how this works, note that if you use the series representation of $$f(x)$$ and square it, you will obtain the square of each of the original terms, but also cross-terms,

$$
f^2(x) = a_0^2 + a_1^2 \cos(\pi x/L) + a_2^2 \cos(2 \pi x/L) + 2 a_0 a_1 \cos(\pi x/L) + a_0 a_2 \cos(2 \pi x/L) + 2 \cos(\pi x/L)  \cos(2 \pi x/L).
$$

However, if we integrated over $$x$$, the cross-terms vanish thanks to the orthogonality relationship. 

If we call the quantity

$$
\frac{1}{L} \int_0^L f^2(x) \, dx
$$

*variance* for now (you might call it something else in a different context), then what Plancherel's theorem says is that we don't need to compute all the cross-terms to compute variance. In fact, each term in the cosine series makes a *unique* contribute to the total variance---this property is critical to ascribing physical meaning.

### Physical interpretation

Consider now a simple equation for a fluid with velocity $$u$$,

$$
u_t - c u_x = 0
$$

that describes waves in a periodic domain of length L. If we multiply the equation by $$u/(2L)$$ and then integrate in space you have that,

$$
\partial_t \frac{1}{2L} \int u^2 dx - \frac{c}{2L} \int u u_x dx = 0.
$$

The last term goes to zero for a periodic domain and thus the quantity
 
$$
\textrm{KE} \equiv \frac{1}{2L} \int u^2 dx
$$

is conserved in time. Assuming that $$u$$ is a velocity, we would typically identify this as the horizontally-averaged kinetic energy. Note that this is essentially the same quantity that we called variance above.

Actual solutions to this wave equation are

$$
u(x,t) = u_n \cos( (n \pi/L) (x+ct) + \phi )
$$

for any $$n=0,1,2,..$$ where $$\phi$$ is the phase of the wave solution at time $$t=0$$. If you take the velocity field $$u(x,t)$$ and project it onto the Fourier modes at time $$t=0$$ (or any time really)

Project onto the Fourier modes, and then you can interpret each mode as representing part of the kinetic energy of the fluid. Our previous interpretation of Plancherel's theorem has even more meaning now because each mode represents a *unique* contribution to the total kinetic energy of the fluid. This uniqueness means that it is okay to ascribe physical meaning to each mode, because each mode can exist in isolation.

It is help to consider an example of where this does *not* hold true.

## Solution Series

Rather than project the equations of motion onto the Fourier series, we can project the equations of motion onto the solutions of the equations. Mathematically this is a pretty minor distinction, but physically,

Each mode in the system 

## Rotating Bousinesq equations in arbitrary stratification


