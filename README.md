GLOceanKit
===========
GLOceanKit is a collection of models and analysis tools for physical oceanography.

The code is in written in two different languages: Matlab and Objective-C, but not all models or analysis tools are available in both languages.

[Matlab](Matlab/)
-------
The Matlab directory contains the following subdirectories of models and tools,
- [Advection-Diffusion Estimation](Matlab/AdvectionDiffusionEstimation) Tools for computing estimating velocity field parameters (strain, vorticity, divergence). From Oscroft, Sykulski & Early (2021).
- [Advection-Diffusion Models](Matlab/AdvectionDiffusionModels) Code for generating particles in advection diffusion models with boundaries.
- [Boussinesq2D](Matlab/Boussinesq2D) 2D nonlinear spectral Boussinesq model.
- [Diffusivity](Matlab/Diffusivity) A collection of analysis tools for computing relative diffusivity from particles.
- [InternalModes](Matlab/InternalModes) Tools solving the vertical mode eigenvalue problem with very high accuracy. From Early, Lelong & Smith (2020).
- [InternalWaveModel](Matlab/InternalWaveModel) A linear internal wave model.
- [InternalWaveSpectrum](Matlab/InternalWaveSpectrum) Tools for computing the Garrett-Munk spectrum and its approximations.
- [OceanBoundaryLayer](Matlab/OceanBoundaryLayer) A few simple ocean boundary layer models taken from Elipot and Gille (2009).
- [Quasigeostrophy](Matlab/Quasigeostrophy) Tools for analyzing the output of the Quasigeostrophic model.

[Objective-C](GLOceanKit/)
-------
Contains internal modes routines, internal wave model, and a QG model.

git-lfs
--------
This repo links to the lfs for some precomputed internal wave modes, but does not download them by default. To override these settings, see [this comment](https://github.com/git-lfs/git-lfs/issues/2717). I think, that if you just do,
`git config lfs.fetchexclude ""`
then it'll remove the exclusion and you can started to download those big files.
