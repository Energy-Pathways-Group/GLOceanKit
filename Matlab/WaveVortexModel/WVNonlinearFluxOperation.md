
A WVNonlinearFluxOperation defines how energy moves between the wave-vortex coefficients---these are the nonlinear terms in the equations of motion, transformed into wave-vortex space. The most basic implementation is a freely evolving, unforced (and undamped) flux. Subclasses of WVNonlinearFluxOperation can implement custom forcing and damping.

The nonlinear terms in the equations of motion are,
$$
\begin{subequations}
\begin{align}
    \textrm{uNL}=& u \partial_x u + v \partial_y u + w \partial_z u \\
    \textrm{vNL}=&u \partial_x v + v \partial_y v + w \partial_z v \\
    \textrm{nNL}=& u \partial_x \eta + v \partial_y \eta + w \left(\partial_z \eta +\eta \partial_z \ln N^2 \right)
\end{align}
\end{subequations}
$$
which, after a transformation into wave-vortex space, define the flux coefficients,
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

The output variables *must* be at least one of {Fp,Fm,F0}, in that order. The properties `doesFluxAp` etc. should be appropriately set to match the output.

You may also optionally output additional variables that are computed as a by product of your flux calculation. Those variables will then be cached, making other function calls much quicker.

- Declaration: classdef WVNonlinearFluxOperation < WVOperation
