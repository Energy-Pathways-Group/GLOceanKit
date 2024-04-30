%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% WaveVortexDecompositionTest
%
% This script tests the API decompose an existing (u,v,w,rho_prime) into
% wave-vortex components
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% April 12th, 2018      Version 1.0

% wvt = WVTransformHydrostatic([15e3, 15e3, 5000], [64 64 33], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
wvt = WVTransformBoussinesq([15e3, 15e3, 5000], [64 64 33], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
% wvt = WVTransformConstantStratification([15e3, 15e3, 5000], [64 64 33]);
wvt.initWithRandomFlow();

% wvt.removeEnergyFromAliasedModes();
% spatialFlux = WVNonlinearFluxSpatial(wvt);
% standardFlux = WVNonlinearFlux(wvt,shouldAntialias=0);
% 
% [SpatialFp,SpatialFm,SpatialF0] = spatialFlux.compute(wvt);
% [StandardFp,StandardFm,StandardF0] = standardFlux.compute(wvt);

% return

%%
[Fp,Fm,F0] = wvt.nonlinearFlux();
% profile on
tic
for i=1:50
    wvt.t = i; % prevent caching
    [Fp,Fm,F0] = wvt.nonlinearFlux();
end
toc

%%
wvt.nonlinearFluxOperation = WVNonlinearFluxSpatial(wvt);
[Fp,Fm,F0] = wvt.nonlinearFlux();

tic
for i=1:50
    wvt.t = i;
    [Fp,Fm,F0] = wvt.nonlinearFlux();
end
toc

% 18:15 2024-02-04 pre-omptimize: Elapsed time is 13.436532 seconds.
% 18:23 2024-02-04 remove squeeze: Elapsed time is 9.491209 s seconds.
% 18:38 2024-02-04 precompute double MM: Elapsed time is 8.681165 seconds.
% 08:09 2024-02-05 loop over kUnique: Elapsed time is 5.783941 seconds.
% 10:00 2024-04-27 [Nj Nkl] and dealiasing refactor: Elapsed time is 2.961629 seconds.

%%
profile on
for i=1:10
    wvt.t = i;
    [Fp,Fm,F0] = wvt.nonlinearFlux();
end
profile viewer
%%

% 18:15 2024-02-04 pre-omptimize: Elapsed time is 3.698 s seconds.
% 18:23 2024-02-04 remove squeeze: Elapsed time is 2.725 s seconds.
% 18:38 2024-02-04 precompute double MM: Elapsed time is 2.289  seconds.
% 08:09 2024-02-05 loop over kUnique: Elapsed time is 1.649 seconds.
% 13:23 2024-04-27 [Nj Nkl] and dealiasing refactor: Elapsed time is 0.886 seconds.
% profile viewer

%%
A = randn(wvt.spectralMatrixSize) + sqrt(-1)*randn(wvt.spectralMatrixSize);
A(wvt.Kh==0) = 0;

profile on
for i=1:100
    % u = wvt.transformToSpatialDomainWithF(A0=A,Apm=A);
    % w = wvt.transformToSpatialDomainWithG(A0=A,Apm=A);
    w = wvt.transformToSpatialDomainWithGUnrolled(A,A);
end
profile viewer