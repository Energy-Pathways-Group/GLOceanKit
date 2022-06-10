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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 128;
aspectRatio = 1;

Lx = 100e3;
Ly = aspectRatio*100e3;
Lz = 5000;

Nx = N;
Ny = aspectRatio*N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1 == 0
    wvm = WaveVortexModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
else
    rho0 = 1025; g = 9.81;
    rho = @(z) -(N0*N0*rho0/g)*z + rho0;
    N2Function = @(z) N0*N0*ones(size(z));
    dLnN2Function = @(z) zeros(size(z));
    wvm = WaveVortexModelHydrostatic([Lx, Ly, Lz], [Nx, Ny, Nz-1], latitude, rho,'N2func', N2Function, 'dLnN2func',dLnN2Function);
end
[ApIO,AmIO,ApIGW,AmIGW,A0G,A0G0,A0rhobar] = wvm.GenerateRandomFlowState();
Ap = ApIO + ApIGW;
Am = AmIO + AmIGW;
A0 = A0G + A0G0 + A0rhobar;

Ubar = wvm.UAp.*Ap + wvm.UAm.*Am + wvm.UA0.*A0;
% Nbar = boussinesq.NAp.*Ap + boussinesq.NAm.*Am + boussinesq.NA0.*A0;
% 
% profile on
% for i=1:100
%     u = boussinesq.transformToSpatialDomainWithF(Ubar);
%     ubar = boussinesq.transformFromSpatialDomainWithF(u);
%     eta = boussinesq.transformToSpatialDomainWithG(Nbar);
%     nbar = boussinesq.transformFromSpatialDomainWithG(eta);
% end
% profile viewer

% profile on
tic
for i=1:15
[Fp,Fm,F0] = wvm.NonlinearFluxAtTimeNoMask(10,Ap,Am,A0);
end
toc

tic
for i=1:15
[Fp,Fm,F0] = wvm.NonlinearFluxAtTime(10,Ap,Am,A0);
end
toc


% profile viewer

