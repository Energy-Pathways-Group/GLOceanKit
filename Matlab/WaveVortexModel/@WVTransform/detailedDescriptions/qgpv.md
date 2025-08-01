- topic: Potential Vorticity & Enstrophy

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
\textrm{QGPV} = \mathcal{DFT}_x^{-1} \left[\mathcal{DFT}_y^{-1} \left[ \mathcal{F}^{-1} \left[ \textrm{QGPV}^{klj} A_0^{klj} \right] \right] \right].
$$

where [$$\textrm{QGPV}^{klj}$$ are linear coefficients](./a0_qgpv_factor.html).
