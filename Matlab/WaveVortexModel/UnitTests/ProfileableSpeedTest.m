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
% wvt = WVTransformConstantStratification([15e3, 15e3, 5000], [128, 128, 65]);
wvt.initWithRandomFlow();

[Fp,Fm,F0] = wvt.nonlinearFlux();
% profile on
tic
for i=1:50
    wvt.t = i; % prevent caching
    [Fp,Fm,F0] = wvt.nonlinearFlux();
end
toc

wvt.nonlinearFluxOperation = WVNonlinearFluxSpatial(wvt);
[Fp,Fm,F0] = wvt.nonlinearFlux();
tic
for i=1:50
    wvt.t = i;
    [Fp,Fm,F0] = wvt.nonlinearFlux();
end
toc

% profile viewer

