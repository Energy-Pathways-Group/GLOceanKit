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

% wvt = WVTransformHydrostatic([15e3, 15e3, 5000], [4, 8, 5], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
% wvt = WVTransformBoussinesq([15e3, 15e3, 5000], [4, 8, 5], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
wvt = WVTransformConstantStratification([15e3, 15e3, 5000], [4, 8, 5]);
wvt.initWithRandomFlow();

% profile on
tic
for i=1:15
[Fp,Fm,F0] = wvt.nonlinearFlux();
end
toc

% profile viewer

