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

boussinesq = Boussinesq3DConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

shouldExcludeNyquist = 1;
Ap = InternalWaveModel.GenerateHermitianRandomMatrix( size(boussinesq.ApU), shouldExcludeNyquist );
Am = InternalWaveModel.GenerateHermitianRandomMatrix( size(boussinesq.ApU), shouldExcludeNyquist );
A0 = InternalWaveModel.GenerateHermitianRandomMatrix( size(boussinesq.ApU), shouldExcludeNyquist );

% Ubar = boussinesq.UAp.*Ap + boussinesq.UAm.*Am + boussinesq.UA0.*A0;
% Nbar = boussinesq.NAp.*Ap + boussinesq.NAm.*Am + boussinesq.NA0.*A0;
% 
% profile on
% for i=1:100
%     u = boussinesq.TransformToSpatialDomainWithF(Ubar);
%     ubar = boussinesq.TransformFromSpatialDomainWithF(u);
%     eta = boussinesq.TransformToSpatialDomainWithG(Nbar);
%     nbar = boussinesq.TransformFromSpatialDomainWithG(eta);
% end
% profile viewer

profile on
for i=1:10
%     [u,v,w,eta] = boussinesq.VelocityField(Ap,Am,A0);
%     [uNL,vNL,nNL] = boussinesq.NonlinearFluxFromSpatial(u,v,w,eta);
%     [uNL,vNL,nNL] = boussinesq.NonlinearFlux(Ap,Am,A0);
%     [Ap,Am,A0] = boussinesq.Project(uNL,vNL,nNL);
[Fp,Fm,F0] = boussinesq.NonlinearFluxAtTime(Ap,Am,A0,10);
end
profile viewer

