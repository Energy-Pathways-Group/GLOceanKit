- Topic: Operations — Transformations
- Declaration: u = transformToSpatialDomainWithF(options)
- Parameter Apm: (optional) variable with dimensions $$(k,l,j)$$ to be transformed with the wave modes
- Parameter A0: (optional) variable with dimensions $$(k,l,j)$$ to be transformed with the geostrophic modes

This is the component of the [inverse discrete transformation](/mathematical-introduction/transformations.html) $$D^{-1}$$ that projects from the vertical modes $F$, followed by a transformation $$ (k,l) \mapsto (x,y)$$ with a discrete Fourier transform. Mathematically we write,

$$
f(x,y,z) =  \mathcal{DFT}_x^{-1} \left[\mathcal{DFT}_y^{-1} \left[ \mathcal{F}^{-1} \left[ \tilde{f}^{klj} \right] \right] \right] 
$$

The $$F$$ mode projection is applicable to dynamical variables $$u$$, $$v$$, $$p$$.

As noted in [Early, et al. (2021)](https://doi.org/10.1017/jfm.2020.995), the vertical transforms $\mathcal{F}$ and $\mathcal{G}$ require a matrix multiplication and thus have a computational cost of,

$$
N_z^2 N_x N_y
$$  

while the horizontal transforms use an FFT algorithm and therefore scale as,

$$
\frac{5}{2} N_z N_x N_y \log_2 N_x N_y - N_z N_x N_y
$$

Assuming that $$N_x = N_y$$, the total computational cost of the horizontal and vertical transforms are approximately equal when $$10 log_2 N_x = N_z$$ , or $$13 log_2 N_x = N_z$$ for the hydrostatic case. This means that for a horizontal resolution of $$256^2$$ the horizontal transformations dominate the computation time until approximately $$80-100$$ vertical modes are used.
