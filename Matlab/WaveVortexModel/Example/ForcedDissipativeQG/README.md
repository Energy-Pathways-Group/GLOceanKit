Forced Dissipative Quasigeostrophic Turbulence
==============

This example shows how to spin up a forced-dissipative QG turbulence simulation, double its resolution, and add particles that track various scalar values along their path.

1. `Spinup.m` Spin up a turbulent fluid using narrow-band forcing and write it to file.
2. `RestartAndDoubleResolution.m` Restart the simulation from existing output with double the resolution.
3. `RestartAndAddParticles.m` Restart the simulation from existing output and add particles, tracking various fluid properties along their trajectories.
4. `DiagnoseEnergyFluxes.m` Diagnose the energy fluxes (inertial, forcing, and damping) as a function of wavenumber for any of the existing output files.