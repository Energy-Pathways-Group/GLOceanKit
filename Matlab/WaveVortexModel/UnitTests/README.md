WaveVortexModel Unit Tests
==============

The WaveVortexModel uses Matlab [unit testing framework](https://www.mathworks.com/help/matlab/write-unit-tests.html) to test various aspects of the model.

In general the unit tests will initialize a small resolution version of each of the three different transforms (constant, hydrostatic, Boussinesq), and then apply the various tests. Test coverage is fairly small at the moment, but the goal is to have all APIs tested. The models are initialized with constant stratification, because there we have analytical solutions that we can test.

- `TestSpectralDifferentiationXY` tests differentiation of the first, second, third and fourth derivatives in each of the x and y directions at each resolved wavenumber. This test calls `diffX` and `diffY`.
- `TestSpectralDifferentiationZ` tests the first four spectral derivatives in the vertical. This test calls `diffZ`.
- `TestOrthogonalSolutionGroups` test each solution from the four solution groups (geostrophic, mda, igw, io) validating (u,v,w,eta,p) as well as the first derivatives of u and eta. This test calls `initWithUVEta`, `u`, `v`, `w`, `eta`, `p`, `transformToSpatialDomainWithFAllDerivatives` and `transformToSpatialDomainWithGAllDerivatives`.  
