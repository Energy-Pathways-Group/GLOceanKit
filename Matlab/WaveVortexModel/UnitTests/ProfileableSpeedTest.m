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

% profile on
% wvt = WVTransformHydrostatic([15e3, 15e3, 5000], 2*[64 64 33], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
wvt = WVTransformHydrostatic([15e3, 15e3, 5000], [1024 1024 30], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
% wvt = WVTransformBoussinesq([15e3, 15e3, 5000], [64 64 33], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
% wvt = WVTransformConstantStratification([15e3, 15e3, 5000], [128 128 128]);
% wvt = WVTransformSingleMode([2000e3 1000e3], 2*[256 128], h=0.8, latitude=25);
% profile viewer

%%
% profile on
wvt.initWithRandomFlow();
% profile viewer

% wvt.removeEnergyFromAliasedModes();
% spatialFlux = WVNonlinearFluxSpatial(wvt);
% standardFlux = WVNonlinearFlux(wvt,shouldAntialias=0);
% 
% [SpatialFp,SpatialFm,SpatialF0] = spatialFlux.compute(wvt);
% [StandardFp,StandardFm,StandardF0] = standardFlux.compute(wvt);

% return

u = wvt.u;
du = zeros(size(u));
u_bar = wvt.transformFromSpatialDomainWithFourier(u);
%%

profile on
tic
for i=1:50
    % du = wvt.diffX(u);
    % du = wvt.fastTransform.diffXIntoArray(u,du);
    u_bar = wvt.transformFromSpatialDomainWithFourier(u);
    % u_bar = u_bar*wvt.fastTransform.dftXY.scaleFactor;
    % u = wvt.transformToSpatialDomainWithFourier(u_bar);
end
toc
profile viewer

%%
% unsorted is fastest!!! DFT sorted is 2x slower, WV sorted is 10% slower
% tic
% for i=1:500
%     wvt.dftBuffer(wvt.dftPrimaryIndex) = wvt.wvBuffer(wvt.wvPrimaryIndex);
%     wvt.dftBuffer(wvt.dftConjugateIndex) = conj(wvt.wvBuffer(wvt.wvConjugateIndex));
% end
% toc
% 
% %%
% % unsorted is 10% slower than sorting either dimension
% tic
% for i=1:500
%     wvt.wvBuffer = wvt.dftBuffer(wvt.dftPrimaryIndex);
% end
% toc

%%
[Fp,Fm,F0] = wvt.nonlinearFlux();

% return
profile on
tic
for i=1:5
    wvt.t = i; % prevent caching
    [Fp,Fm,F0] = wvt.nonlinearFlux();
end
toc
profile viewer

return

%%
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
% 16:09 2024-05-01 Improved memory reshuffle: Elapsed time is 1.784329 seconds.

%%
profile on
for i=1:200
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
% 2024-04-30 Finally got this down to something reasonable!!!
% 
A = randn(wvt.spectralMatrixSize) + sqrt(-1)*randn(wvt.spectralMatrixSize);
A(wvt.Kh==0) = 0;

% profile on
tic
for i=1:500
    u = wvt.transformToSpatialDomainWithF(A0=A,Apm=A);
    w = wvt.transformToSpatialDomainWithG(A0=A,Apm=A);
    % w = wvt.transformToSpatialDomainWithGUnrolled(A,A);
end
toc
% profile viewer
% 2.7, 2.8

%%
% 2024-04-30 Finally got this down to something reasonable!!!
% 
A = randn(wvt.spatialMatrixSize);

% profile on
tic
for i=1:5000
    u = wvt.transformFromSpatialDomainWithFourier(A);
end
toc
% profile viewer
% 2.9, 2.8

%%
A = randn(100,100);
B = randn(100,1000) + sqrt(-1)*randn(100,1000);
A(1,:)=0; A(:,1)=0;
tic;
for i=1:100
    C = A*B;
end
toc
