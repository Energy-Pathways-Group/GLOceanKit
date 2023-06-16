- topic: State Variables

The quasigeostrophic potential vorticty (QGPV) is defined as,

$$
\textrm{QGPV} \equiv \partial_x v - \partial_y u - f_0 \partial_z \eta_\textrm{e}
$$

where $$\eta_\textrm{e}$$ is the linear approximation to the isopycnal displacement, and related to excess density with $$N^2 \eta_\textrm{e}= \frac{g}{\rho_0} \rho$$.

This same quantity can be computed from a stream function $$\psi$$ with,

$$
\textrm{QGPV} = \nabla^2 \psi + \frac{d}{dz}\left( \frac{f^2}{N^2} \frac{d \psi}{dz} \right)
$$

The coefficients $$A_0$$ are linearly related to the QGPV and the streamfunction such that,

$$
\textrm{QGPV} = \mathcal{DFT}_x^{-1} \left[\mathcal{DFT}_y^{-1} \left[ \mathcal{F}^{-1} \left[-\left( K^2 + \frac{f^2}{gh_j} \right)\frac{g}{f} A_0^{klj} \right] \right] \right].
$$

For hydrostatic motions, the coefficient can be related to the frequency of the internal gravity waves,

$$
\frac{\omega^2}{h f_0} = \left( K^2 + \frac{f^2}{gh} \right)\frac{g}{f}
$$

and thus, in practice, the implementation in the hydrostatic transform uses,

$$
\textrm{QGPV} = \mathcal{DFT}_x^{-1} \left[\mathcal{DFT}_y^{-1} \left[ \mathcal{F}^{-1} \left[-\frac{\omega_j^2}{h f_0} A_0^{klj} \right] \right] \right].
$$
